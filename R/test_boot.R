#' Bootstrap median S_norm (percentile and pivotal CIs)
#'
#' Resamples subjects with replacement, refits the selected NLME model on
#' control + monotherapy data, recomputes subject-wise signed standardized
#' distances \eqn{S_i}, and summarizes the bootstrap distribution of the
#' median \eqn{S_i}.
#'
#' \deqn{S_i = \mathrm{sign}\!\left(\overline{y_0 - y}\right) \cdot D_i / n_i}
#'
#' @param fit synpdx_fit
#' @param data data.frame with columns subject,time,treatment,logvol_norm
#' @param control character. Control arm name.
#' @param drug_a character. Drug A arm name.
#' @param drug_b character. Drug B arm name.
#' @param combo  character. Combination arm name.
#' @param B integer. Number of bootstrap replicates.
#' @param seed integer. RNG seed.
#' @param saem_boot list. Overrides for \code{saemix::saemixControl()} in bootstrap refits.
#' @param out_dir character or NULL. Optional output directory for CSVs and PNG.
#'
#' @return list with elements:
#'   \item{summary}{data.frame with percentile and pivotal 95% CIs and interaction labels}
#'   \item{subject_Si}{data.frame of per-replicate subject \eqn{S_i}}
#'   \item{plot}{ggplot showing bootstrap density and CI bands}
#' @export

synpdx_bootstrap <- function(fit, data, control, drug_a, drug_b, combo,
                             B = 100, seed = 12345, saem_boot = list(), out_dir = NULL) {
  stopifnot(inherits(fit, "synpdx_fit"))
  set.seed(seed)

  # Observed S_norm via compute_Di
  obs <- synpdx_compute_Di(fit, data, control, drug_a, drug_b, combo, cutoff = "auto")
  combo_df <- obs$combo_df
  Di_vec   <- obs$Di_vec

  combo_n_obs_original <- combo_df |>
    dplyr::group_by(.data$subject) |>
    dplyr::summarise(n_obs = dplyr::n(), .groups = "drop")

  result_tbl_obs <- combo_df |>
    dplyr::mutate(diff = .data$predicted_logvol - .data$DV) |>
    dplyr::group_by(.data$subject) |>
    dplyr::summarise(mean_resid = mean(.data$diff, na.rm = TRUE), .groups = "drop") |>
    dplyr::left_join(tibble::tibble(subject = names(Di_vec), Di = as.numeric(Di_vec)),
                     by = "subject") |>
    dplyr::left_join(combo_n_obs_original, by = "subject") |>
    dplyr::filter(is.finite(.data$Di), is.finite(.data$n_obs), .data$n_obs > 0) |>
    dplyr::mutate(S_norm = sign(.data$mean_resid) * (.data$Di / .data$n_obs))

  S_norm_median_original <- stats::median(result_tbl_obs$S_norm, na.rm = TRUE)

  data_processed <- fit$data_processed
  arms_all <- c(control, drug_a, drug_b, combo)
  subjects_all <- unique(as.character(data_processed$subject))

  selected_model <- fit$selected_model
  model_fun <- selected_model@model@model

  # Helpers (all qualified to avoid missing-import NOTES)
  helper_make_row <- function(x, n) {
    if (n == 0) return(matrix(numeric(0), nrow = 0, ncol = length(x)))
    xv <- as.numeric(x); cn <- names(x)
    M  <- matrix(rep(xv, times = n), nrow = n, byrow = TRUE)
    if (!is.null(cn) && length(cn) == ncol(M)) colnames(M) <- cn
    M
  }
  pred_arm <- function(times, IA, IB, logpar_row){
    n <- length(times); if (n == 0) return(numeric(0))
    psi <- helper_make_row(logpar_row, n)
    as.numeric(model_fun(psi = psi, id = rep(1L, n), xidep = data.frame(time = times, I_A = IA, I_B = IB)))
  }
  null_pred <- function(times, logpar_row){
    if (length(times) == 0) return(numeric(0))
    yA <- pred_arm(times, 1L, 0L, logpar_row)
    yB <- pred_arm(times, 0L, 1L, logpar_row)
    yC <- pred_arm(times, 0L, 0L, logpar_row)
    yA + yB - yC
  }
  null_jac <- function(times, logpar_row, eta_names, h_scale = 1e-3){
    n <- length(times); q <- length(eta_names)
    if (q == 0 || n == 0) return(matrix(0, n, 0))
    J <- matrix(0, n, q, dimnames = list(NULL, eta_names))
    base <- logpar_row
    for (k in seq_len(q)) {
      nm <- eta_names[k]
      step <- h_scale * max(1, abs(base[[nm]]))
      p1 <- base; p1[[nm]] <- base[[nm]] + step
      m1 <- base; m1[[nm]] <- base[[nm]] - step
      y_plus  <- null_pred(times, p1)
      y_minus <- null_pred(times, m1)
      J[, k] <- (y_plus - y_minus) / (2 * step)
    }
    J
  }
  quad_form <- function(C, d) {
    val <- try({ L <- chol(C); z <- backsolve(L, forwardsolve(t(L), d)); sum(z^2) }, silent=TRUE)
    if (!inherits(val,"try-error") && is.finite(val)) return(as.numeric(val))
    x <- try(qr.solve(C, d), silent=TRUE)
    if (!inherits(x,"try-error") && all(is.finite(x))) return(as.numeric(crossprod(d, x)))
    pinv <- try(MASS::ginv(C), silent=TRUE)
    if (!inherits(pinv,"try-error") && all(is.finite(pinv))) return(as.numeric(t(d) %*% pinv %*% d))
    NA_real_
  }

  boot_medians <- numeric(0)
  all_subject_Si <- list()

  boot_ctrl <- do.call(saemix::saemixControl, utils::modifyList(
    list(print=FALSE, map=TRUE, fim=FALSE, ll.is=FALSE, nbiter.saemix = c(400,150), nb.chains=3, seed=13579),
    saem_boot
  ))

  for (b in seq_len(B)) {
    # Resample subjects
    boot_subjects <- sample(subjects_all, length(subjects_all), replace = TRUE)
    boot_data <- dplyr::bind_rows(lapply(seq_along(boot_subjects), function(i) {
      s <- boot_subjects[i]
      dplyr::mutate(dplyr::filter(data_processed, .data$subject == s),
                    subject = paste0("boot", b, "_", i))
    })) |> dplyr::mutate(subject = as.character(.data$subject))

    # Refit on control + monotherapies
    boot_fit_data <- dplyr::filter(boot_data, .data$treatment != combo)
    boot_sdm <- saemix::saemixData(
      name.data = boot_fit_data, name.group = "subject",
      name.predictors = c("time","I_A","I_B"), name.response = "logvol_norm"
    )
    fit_try <- try(saemix::saemix(selected_model@model, boot_sdm, control = boot_ctrl), silent = TRUE)
    if (inherits(fit_try, "try-error")) next
    fit_boot <- fit_try

    # Subject-specific cutoff
    cutoff_boot <- boot_data |>
      dplyr::filter(.data$treatment %in% arms_all) |>
      dplyr::group_by(.data$subject, .data$treatment) |>
      dplyr::summarise(tmax_arm = max(.data$time), .groups = "drop") |>
      dplyr::group_by(.data$subject) |>
      dplyr::summarise(cutoff = min(.data$tmax_arm), .groups = "drop")

    # Combo observations up to cutoff
    combo_df_boot <- boot_data |>
      dplyr::filter(.data$treatment == combo) |>
      dplyr::rename(DV = .data$logvol_norm) |>
      dplyr::select(.data$subject, .data$time, .data$DV) |>
      dplyr::left_join(cutoff_boot, by = "subject") |>
      dplyr::filter(.data$time <= .data$cutoff) |>
      dplyr::arrange(.data$subject, .data$time)
    if (nrow(combo_df_boot) == 0) next

    # Extract parameters and structures
    par_model <- as.character(fit_boot@model@name.modpar)
    mu_full   <- fit_boot@results@fixed.effects
    if (is.null(names(mu_full)) || any(names(mu_full) == "")) names(mu_full) <- par_model[seq_len(length(mu_full))]
    Omega_full <- fit_boot@results@omega
    if (is.null(rownames(Omega_full)) || is.null(colnames(Omega_full))) {
      rn <- par_model[seq_len(nrow(Omega_full))]; rownames(Omega_full) <- colnames(Omega_full) <- rn
    }
    par_common <- Reduce(intersect, list(par_model, names(mu_full), colnames(Omega_full)))
    if (length(par_common) == 0L) next
    mu_vec  <- as.numeric(mu_full[par_common]); names(mu_vec) <- par_common
    Omega   <- Omega_full[par_common, par_common, drop = FALSE]
    eta_names <- names(which(diag(Omega) > .Machine$double.eps))
    Omega_use <- if (length(eta_names)>0) Omega[eta_names, eta_names, drop = FALSE] else matrix(0,0,0)

    ebes <- fit_boot@results@map.psi
    if (is.null(rownames(ebes)) || any(rownames(ebes) == "")) {
      phiM <- fit_boot@results@phiM
      if (!is.null(phiM$group) && length(phiM$group) == nrow(ebes))
        rownames(ebes) <- as.character(phiM$group)
    }
    ebes <- ebes[, intersect(colnames(ebes), names(mu_vec)), drop = FALSE]

    respar <- as.numeric(fit_boot@results@respar)
    sigma2_model <- if (length(respar) >= 1L) respar[1]^2 else NA_real_

    fit_df <- boot_fit_data |>
      dplyr::transmute(subject = as.character(.data$subject),
                       time = .data$time,
                       arm  = .data$treatment,
                       DV   = .data$logvol_norm) |>
      dplyr::left_join(cutoff_boot, by = "subject") |>
      dplyr::filter(.data$time <= .data$cutoff)

    # Empirical residual variance on refit
    res_all <- c()
    for (s in unique(fit_df$subject)) {
      mu_i <- mu_vec
      if (s %in% rownames(ebes)) {
        eta_i <- ebes[s, names(mu_vec), drop = TRUE]; eta_i[is.na(eta_i)] <- 0
        mu_i <- mu_vec + eta_i
      }
      dfs <- dplyr::filter(fit_df, .data$subject == s)
      if (nrow(dfs) == 0) next
      dC <- dplyr::arrange(dplyr::filter(dfs, .data$arm == control), .data$time)
      dA <- dplyr::arrange(dplyr::filter(dfs, .data$arm == drug_a),  .data$time)
      dB <- dplyr::arrange(dplyr::filter(dfs, .data$arm == drug_b),  .data$time)
      if (nrow(dC) > 0) res_all <- c(res_all, dC$DV - pred_arm(dC$time, 0L, 0L, mu_i))
      if (nrow(dA) > 0) res_all <- c(res_all, dA$DV - pred_arm(dA$time, 1L, 0L, mu_i))
      if (nrow(dB) > 0) res_all <- c(res_all, dB$DV - pred_arm(dB$time, 0L, 1L, mu_i))
    }
    sigma2_emp  <- stats::var(res_all, na.rm = TRUE)
    sigma2_eff <- if (is.finite(sigma2_emp) && sigma2_emp > 0) sigma2_emp else sigma2_model
    if (!is.finite(sigma2_eff) || sigma2_eff <= 0) next

    # Subject-wise S_i and median
    S_vec <- c(); subj_rows <- list()
    for (sid in unique(combo_df_boot$subject)) {
      dsub <- dplyr::arrange(dplyr::filter(combo_df_boot, .data$subject == sid), .data$time)
      y_obs <- dsub$DV; times <- dsub$time; n_i <- length(times)
      if (any(!is.finite(y_obs)) || n_i < 2) next

      mu_i <- mu_vec
      if (sid %in% rownames(ebes)) {
        eta_i <- ebes[sid, names(mu_vec), drop = TRUE]; eta_i[is.na(eta_i)] <- 0
        mu_i <- mu_vec + eta_i
      }
      J <- null_jac(times, mu_i, eta_names)
      y0 <- null_pred(times, mu_i)
      d  <- y_obs - y0
      Cov_i <- if (length(eta_names) > 0) J %*% Omega_use %*% t(J) + diag(sigma2_eff, n_i) else diag(sigma2_eff, n_i)
      Di_val <- quad_form(Cov_i, d); if (!is.finite(Di_val)) next
      mean_resid <- mean(y0 - y_obs, na.rm = TRUE)
      S_i <- sign(mean_resid) * (Di_val / n_i)
      S_vec <- c(S_vec, S_i)
      subj_rows[[length(subj_rows)+1]] <- tibble::tibble(
        replicate = b, subject = sid, n_obs = n_i,
        mean_resid = mean_resid, Di = Di_val, S_norm = S_i
      )
    }
    if (length(S_vec) > 0 && any(is.finite(S_vec))) {
      boot_medians <- c(boot_medians, stats::median(S_vec, na.rm = TRUE))
      all_subject_Si[[length(all_subject_Si)+1]] <- dplyr::bind_rows(subj_rows)
    }
  }

  # CIs and interaction labels
  if (length(boot_medians) > 0) {
    qs <- stats::quantile(boot_medians, c(0.025, 0.975), type = 7, na.rm = TRUE)
    qL <- as.numeric(qs[1]); qU <- as.numeric(qs[2])

    # Pivotal: 2*theta_hat - percentile
    ci_lower_piv <- 2*S_norm_median_original - qU
    ci_upper_piv <- 2*S_norm_median_original - qL
    interaction_piv <- if (ci_lower_piv > 0) "synergistic" else if (ci_upper_piv < 0) "antagonistic" else "additive"

    # Percentile: [qL, qU]
    ci_lower_perc <- qL
    ci_upper_perc <- qU
    interaction_perc <- if (ci_lower_perc > 0) "synergistic" else if (ci_upper_perc < 0) "antagonistic" else "additive"
  } else {
    ci_lower_piv <- ci_upper_piv <- ci_lower_perc <- ci_upper_perc <- NA_real_
    interaction_piv <- interaction_perc <- "Failed"
  }

  # Plot
  plot_df <- tibble::tibble(boot_medians = boot_medians)
  dobj <- if (nrow(plot_df) > 1 && all(is.finite(plot_df$boot_medians))) {
    dens <- stats::density(plot_df$boot_medians, na.rm = TRUE)
    tibble::tibble(x = dens$x, y = dens$y)
  } else NULL
  band_height <- if (!is.null(dobj)) max(dobj$y) * 0.15 else 1

  p <- ggplot2::ggplot() +
    { if (!is.null(dobj)) ggplot2::geom_area(data = dobj, ggplot2::aes(x = .data$x, y = .data$y), alpha = 0.15) } +
    ggplot2::geom_histogram(
      data = plot_df,
      ggplot2::aes(x = .data$boot_medians, y = ggplot2::after_stat(density)),
      bins = 20, fill = NA, color = "black"
    ) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = 0, linetype = "Null (0)"), linewidth = 1) +
    { if (is.finite(ci_lower_perc) && is.finite(ci_upper_perc))
      ggplot2::annotate("rect", xmin = ci_lower_perc, xmax = ci_upper_perc,
                        ymin = 0, ymax = band_height, alpha = 0.35) } +
    { if (is.finite(ci_lower_piv) && is.finite(ci_upper_piv))
      ggplot2::annotate("rect", xmin = ci_lower_piv, xmax = ci_upper_piv,
                        ymin = 0, ymax = band_height, alpha = 0.35) } +
    ggplot2::scale_linetype_manual(name = NULL, breaks = "Null (0)", values = "dashed") +
    ggplot2::guides(linetype = ggplot2::guide_legend(override.aes = list(linewidth = 1))) +
    ggplot2::labs(title = "Bootstrap distribution of median S_norm",
                  subtitle = "Shaded: Percentile and Pivotal 95% CIs",
                  x = "Median of S_i / n_i", y = "Density") +
    ggplot2::theme_bw(base_size = 18) +
    ggplot2::theme(legend.position = "bottom",
                   panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1.2))

  res_tbl <- dplyr::bind_rows(
    tibble::tibble(method="Pivotal Bootstrap",   n_boot=length(boot_medians),
                   S_norm_median_observed=S_norm_median_original,
                   CI_95_lower=ci_lower_piv, CI_95_upper=ci_upper_piv,
                   overall_interaction=interaction_piv),
    tibble::tibble(method="Percentile Bootstrap", n_boot=length(boot_medians),
                   S_norm_median_observed=S_norm_median_original,
                   CI_95_lower=ci_lower_perc, CI_95_upper=ci_upper_perc,
                   overall_interaction=interaction_perc)
  )
  subj_Si_df <- if (length(all_subject_Si) > 0) dplyr::bind_rows(all_subject_Si) else tibble::tibble()

  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    utils::write.csv(res_tbl, file.path(out_dir, "bootstrap_summary_result.csv"), row.names = FALSE)
    if (nrow(subj_Si_df)) utils::write.csv(subj_Si_df, file.path(out_dir, "bootstrap_subject_Si.csv"), row.names = FALSE)
    ggplot2::ggsave(file.path(out_dir, "bootstrap_ci_overlay.png"), plot = p, width = 9, height = 7, dpi = 300)
  }
  list(summary = res_tbl, subject_Si = subj_Si_df, plot = p)
}

