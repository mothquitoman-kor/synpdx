# R/test_boot.R
#' Bootstrap median S_norm (percentile and pivotal CIs)
#' @export
synpdx_bootstrap <- function(fit, data, control, drug_a, drug_b, combo,
                             B = 100, seed = 202, cutoff = c("auto","none"),
                             out_dir = NULL) {
  cutoff <- match.arg(cutoff)
  stopifnot(inherits(fit, "synpdx_fit"))
  data <- dplyr::mutate(data, subject = as.character(.data$subject))

    model_name <- fit$selected_model_name
  if (is.null(model_name) || !nzchar(model_name)) {
    stop("`fit$selected_model_name` not found. Pass a fit created by synpdx_fit(model=...).")
  }

  obs_Di <- synpdx_compute_Di(fit, data, control, drug_a, drug_b, combo,
                              cutoff = cutoff, out_dir = NULL)
  combo_df_obs <- dplyr::mutate(obs_Di$combo_df,
                                diff = .data$predicted_logvol - .data$DV)
  n_obs_tbl <- combo_df_obs |>
    dplyr::group_by(.data$subject) |>
    dplyr::summarise(n_i = dplyr::n(),
                     mean_resid = mean(.data$diff, na.rm = TRUE), .groups = "drop")
  Di_df <- tibble::tibble(subject = names(obs_Di$Di_vec),
                          Di = as.numeric(obs_Di$Di_vec))
  obs_tbl <- dplyr::inner_join(n_obs_tbl, Di_df, by = "subject") |>
    dplyr::mutate(S_norm = sign(.data$mean_resid) * .data$Di / .data$n_i)
  S_norm_median_original <- stats::median(obs_tbl$S_norm, na.rm = TRUE)

  subjects_all <- sort(unique(data$subject))
  n_subj <- length(subjects_all); if (n_subj < 2L) stop("Bootstrap needs >= 2 subjects.")
  set.seed(seed)
  boot_medians <- numeric(B)
  all_subject_Si <- vector("list", B)

  for (b in seq_len(B)) {
    cat(sprintf("[bootstrap %s] %d/%d\n", model_name, b, B)); utils::flush.console()

    draws <- sample(subjects_all, size = n_subj, replace = TRUE)
    boot_df <- lapply(seq_along(draws), function(k) {
      s <- draws[k]
      df_s <- dplyr::filter(data, .data$subject == s)
      df_s$subject <- paste0(s, "_boot", b, "_", k)
      df_s
    }) |> dplyr::bind_rows()

    fit_b <- synpdx_model_fit(boot_df, control, drug_a, drug_b, combo,
                              seed = seed + b, model = model_name)

    Di_b  <- synpdx_compute_Di(fit_b, boot_df, control, drug_a, drug_b, combo,
                               cutoff = cutoff, out_dir = NULL)

    combo_df_b <- dplyr::mutate(Di_b$combo_df,
                                diff = .data$predicted_logvol - .data$DV)
    n_tbl_b <- combo_df_b |>
      dplyr::group_by(.data$subject) |>
      dplyr::summarise(n_i = dplyr::n(),
                       mean_resid = mean(.data$diff, na.rm = TRUE), .groups = "drop")
    Di_df_b <- tibble::tibble(subject = names(Di_b$Di_vec),
                              Di = as.numeric(Di_b$Di_vec))
    res_b <- dplyr::inner_join(n_tbl_b, Di_df_b, by = "subject") |>
      dplyr::mutate(S_norm = sign(.data$mean_resid) * .data$Di / .data$n_i,
                    replicate = b, model = model_name) |>
      dplyr::arrange(.data$subject)

    boot_medians[b] <- stats::median(res_b$S_norm, na.rm = TRUE)
    all_subject_Si[[b]] <- res_b[, c("replicate","model","subject","n_i","Di","S_norm")]
  }

  ci_lower_perc <- stats::quantile(boot_medians, 0.025, names = FALSE, type = 7, na.rm = TRUE)
  ci_upper_perc <- stats::quantile(boot_medians, 0.975, names = FALSE, type = 7, na.rm = TRUE)
  diffs <- boot_medians - S_norm_median_original
  ci_lower_piv <- S_norm_median_original - stats::quantile(diffs, 0.975, names = FALSE, type = 7, na.rm = TRUE)
  ci_upper_piv <- S_norm_median_original - stats::quantile(diffs, 0.025, names = FALSE, type = 7, na.rm = TRUE)
  interaction_perc <- if (ci_lower_perc > 0) "synergistic" else if (ci_upper_perc < 0) "antagonistic" else "additive"
  interaction_piv  <- if (ci_lower_piv  > 0) "synergistic" else if (ci_upper_piv  < 0) "antagonistic" else "additive"

  dens_df <- tibble::tibble(x = boot_medians)
  p <- ggplot2::ggplot(dens_df, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_density(adjust = 1) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = S_norm_median_original, linetype = "dashed") +
    ggplot2::labs(x = "Bootstrap median S_norm", y = "Density")

  res_tbl <- tibble::tibble(
    method = c("percentile","pivotal"),
    model  = model_name,
    S_norm_median = S_norm_median_original,
    CI_95_lower = c(ci_lower_perc, ci_lower_piv),
    CI_95_upper = c(ci_upper_perc, ci_upper_piv),
    overall_interaction = c(interaction_perc, interaction_piv)
  )
  subj_Si_df <- if (length(all_subject_Si)) dplyr::bind_rows(all_subject_Si) else tibble::tibble()

  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    utils::write.csv(res_tbl,   file.path(out_dir, "bootstrap_summary_result.csv"), row.names = FALSE)
    if (nrow(subj_Si_df))
      utils::write.csv(subj_Si_df, file.path(out_dir, "bootstrap_subject_Si.csv"), row.names = FALSE)
    ggplot2::ggsave(file.path(out_dir, "bootstrap_ci_overlay.png"), plot = p, width = 9, height = 7, dpi = 300)
  }
  list(summary = res_tbl, subject_Si = subj_Si_df, plot = p)
}
