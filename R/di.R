#' Compute per-subject Mahalanobis distances under Bliss additivity
#' @param fit synpdx_fit
#' @param data data.frame with columns subject,time,treatment,logvol_norm
#' @param control control arm name
#' @param drug_a drug A arm name
#' @param drug_b drug B arm name
#' @param combo combination arm name
#' @param cutoff "auto" or "none"
#' @param out_dir optional directory to write RDS
#' @return list(Di_vec, combo_df, sigma2_used)
#' @export
synpdx_compute_Di <- function(fit, data, control, drug_a, drug_b, combo,
                              cutoff = c("auto","none"), out_dir = NULL) {
  cutoff <- match.arg(cutoff)
  stopifnot(inherits(fit, "synpdx_fit"))
  selected_model <- fit$selected_model
  data_processed <- fit$data_processed

  model_fun <- selected_model@model@model
  par_model <- as.character(selected_model@model@name.modpar)
  mu_full <- selected_model@results@fixed.effects
  if (is.null(names(mu_full)) || any(names(mu_full) == "")) names(mu_full) <- par_model[seq_len(length(mu_full))]
  Omega_full <- selected_model@results@omega
  if (is.null(rownames(Omega_full)) || is.null(colnames(Omega_full))) {
    rn <- par_model[seq_len(nrow(Omega_full))]
    rownames(Omega_full) <- colnames(Omega_full) <- rn
  }
  par_common <- Reduce(intersect, list(par_model, names(mu_full), colnames(Omega_full)))
  mu_vec <- mu_full[par_common]
  Omega  <- Omega_full[par_common, par_common, drop = FALSE]
  eta_names <- names(which(diag(Omega) > .Machine$double.eps))
  Omega_use <- Omega[eta_names, eta_names, drop = FALSE]

  respar <- as.numeric(selected_model@results@respar)
  sigma2_model <- if (length(respar) >= 1L) respar[1]^2 else NA_real_

  make_row <- function(x, n) {
    xv <- as.numeric(x); names(xv) <- names(x)
    matrix(rep(xv, times = n), nrow = n, byrow = TRUE, dimnames = list(NULL, names(xv)))
  }
  pred_arm <- function(times, IA, IB, logpar_row){
    n <- length(times)
    psi <- make_row(logpar_row, n)
    as.numeric(model_fun(psi = psi, id = rep(1L, n),
                         xidep = data.frame(time = times, I_A = IA, I_B = IB)))
  }
  null_pred <- function(times, logpar_row){
    yA <- pred_arm(times, 1L, 0L, logpar_row)
    yB <- pred_arm(times, 0L, 1L, logpar_row)
    yC <- pred_arm(times, 0L, 0L, logpar_row)
    yA + yB - yC
  }
  null_jac_num <- function(times, logpar_row, eta_names, h_scale = 1e-3){
    n <- length(times); q <- length(eta_names)
    J <- matrix(0, n, q, dimnames = list(NULL, eta_names))
    base <- logpar_row
    for (k in seq_len(q)) {
      nm <- eta_names[k]
      step <- h_scale * (1 + abs(base[[nm]]))
      p1 <- base; p1[[nm]] <- base[[nm]] + step
      m1 <- base; m1[[nm]] <- base[[nm]] - step
      y_plus  <- null_pred(times, p1)
      y_minus <- null_pred(times, m1)
      J[, k] <- (y_plus - y_minus) / (2 * step)
    }
    J
  }

  arms_all <- c(control, drug_a, drug_b, combo)
  cutoff_data <- data_processed |>
    dplyr::filter(.data$treatment %in% arms_all) |>
    dplyr::group_by(.data$subject, .data$treatment) |>
    dplyr::summarise(tmax_arm = max(.data$time), .groups = "drop") |>
    dplyr::group_by(.data$subject) |>
    dplyr::summarise(cutoff = min(.data$tmax_arm), .groups = "drop") |>
    dplyr::mutate(subject = as.character(.data$subject))

  combo_df_base <- data_processed |>
    dplyr::filter(.data$treatment == combo) |>
    dplyr::transmute(subject = as.character(.data$subject), time = .data$time, DV = .data$logvol_norm) |>
    dplyr::left_join(cutoff_data, by = "subject") |>
    dplyr::filter(.data$time <= .data$cutoff) |>
    dplyr::arrange(.data$subject, .data$time)

  ebes <- selected_model@results@map.psi
  if (is.null(rownames(ebes)) || any(rownames(ebes) == "")) {
    phiM <- selected_model@results@phiM
    if (!is.null(phiM$group) && length(phiM$group) == nrow(ebes))
      rownames(ebes) <- as.character(phiM$group)
  }
  ebes <- ebes[, intersect(colnames(ebes), names(mu_vec)), drop = FALSE]

  fit_df <- data_processed |>
    dplyr::filter(.data$treatment %in% c(control, drug_a, drug_b)) |>
    dplyr::transmute(subject = as.character(.data$subject), time = .data$time,
                     arm = .data$treatment, DV = .data$logvol_norm) |>
    dplyr::left_join(cutoff_data, by = "subject") |>
    dplyr::filter(.data$time <= .data$cutoff)

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
    if (nrow(dC) > 1) res_all <- c(res_all, dC$DV - pred_arm(dC$time, 0L, 0L, mu_i))
    if (nrow(dA) > 1) res_all <- c(res_all, dA$DV - pred_arm(dA$time, 1L, 0L, mu_i))
    if (nrow(dB) > 1) res_all <- c(res_all, dB$DV - pred_arm(dB$time, 0L, 1L, mu_i))
  }
  sigma2_emp  <- stats::var(res_all, na.rm = TRUE)
  sigma2_y    <- if (is.finite(sigma2_emp) && sigma2_emp > 0) sigma2_emp else sigma2_model

  subjects <- unique(combo_df_base$subject)
  Di_list <- vector("list", length(subjects)); names(Di_list) <- subjects
  for (s in subjects) {
    df <- dplyr::arrange(dplyr::filter(combo_df_base, .data$subject == s), .data$time)
    times <- df$time; y_obs <- df$DV; n_i <- length(times)
    mu_i <- mu_vec
    if (s %in% rownames(ebes)) {
      eta_i <- ebes[s, names(mu_vec), drop = TRUE]; eta_i[is.na(eta_i)] <- 0
      mu_i <- mu_vec + eta_i
    }
    J  <- null_jac_num(times, mu_i, eta_names)
    y0 <- null_pred(times, mu_i)
    d  <- y_obs - y0
    Cov_i <- J %*% Omega_use %*% t(J) + diag(sigma2_y, n_i)
    Di_val <- tryCatch({
      L <- chol(Cov_i); z <- backsolve(L, forwardsolve(t(L), d)); sum(z^2)
    }, error = function(e) {
      eig_min <- suppressWarnings(min(eigen(Cov_i, symmetric = TRUE, only.values = TRUE)$values, na.rm = TRUE))
      add_lambda <- if (is.finite(eig_min)) max(1e-12, 1e-10 * mean(diag(Cov_i)) - eig_min) else 1e-12
      Cov_i2 <- Cov_i + diag(add_lambda, n_i)
      as.numeric(t(d) %*% solve(Cov_i2, d))
    })
    Di_list[[s]] <- Di_val
  }
  Di_vec <- unlist(Di_list)

  pred_additive <- combo_df_base |>
    dplyr::group_by(.data$subject) |>
    dplyr::reframe(
      time   = .data$time,
      DV     = .data$DV,
      cutoff = .data$cutoff,
      predicted_logvol = {
        mu_i <- mu_vec
        if (unique(.data$subject) %in% rownames(ebes)) {
          eta_i <- ebes[unique(.data$subject), names(mu_vec), drop = TRUE]; eta_i[is.na(eta_i)] <- 0
          mu_i <- mu_vec + eta_i
        }
        null_pred(.data$time, mu_i)
      }
    ) |>
    dplyr::ungroup()

  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    saveRDS(Di_vec,        file.path(out_dir, "Di_vec.rds"))
    saveRDS(pred_additive, file.path(out_dir, "combo_df.rds"))
  }
  list(Di_vec = Di_vec, combo_df = pred_additive, sigma2_used = sigma2_y)
}
