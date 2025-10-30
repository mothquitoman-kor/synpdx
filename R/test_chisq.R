# test_chisq.R
#' Group-level chi-square test with AD GOF on PIT
#' @export
synpdx_chisq_test <- function(Di_vec, combo_df, alpha = 0.05, out_dir = NULL) {
  need_cols <- c("subject","time","DV","predicted_logvol","cutoff")
  stopifnot(all(need_cols %in% names(combo_df)))
  stopifnot(!is.null(names(Di_vec)))

  combo_df <- dplyr::mutate(combo_df, subject = as.character(.data$subject)) |>
    dplyr::arrange(.data$subject, .data$time)
  Di_df    <- tibble::tibble(subject = as.character(names(Di_vec)), Di = as.numeric(Di_vec))
  common_subjects <- intersect(unique(combo_df$subject), Di_df$subject)
  stopifnot(length(common_subjects) > 0)
  combo_df <- dplyr::filter(combo_df, .data$subject %in% common_subjects)
  Di_df    <- dplyr::filter(Di_df,    .data$subject %in% common_subjects)

  combo_n_obs <- combo_df |>
    dplyr::group_by(.data$subject) |>
    dplyr::summarise(n_obs = dplyr::n(), .groups = "drop")

  D_total  <- sum(Di_df$Di, na.rm = TRUE)
  df_total <- sum(combo_n_obs$n_obs[match(Di_df$subject, combo_n_obs$subject)], na.rm = TRUE)
  p_value_total <- stats::pchisq(D_total, df = df_total, lower.tail = FALSE)

  result_tbl <- combo_df |>
    dplyr::mutate(diff = .data$predicted_logvol - .data$DV) |>
    dplyr::group_by(.data$subject) |>
    dplyr::summarise(mean_resid = mean(.data$diff, na.rm = TRUE), .groups = "drop") |>
    dplyr::inner_join(Di_df, by = "subject") |>
    dplyr::inner_join(combo_n_obs, by = "subject") |>
    dplyr::mutate(
      sign_resid      = sign(.data$mean_resid),
      Si              = .data$sign_resid * sqrt(.data$Di),
      chisq_threshold = stats::qchisq(0.95, df = .data$n_obs),
      significant     = .data$Di > .data$chisq_threshold,
      interaction     = dplyr::case_when(
        .data$significant & (.data$sign_resid > 0) ~ "synergistic",
        .data$significant & (.data$sign_resid < 0) ~ "antagonistic",
        TRUE                                       ~ "additive"
      ),
      PIT             = stats::pchisq(.data$Di, df = .data$n_obs, lower.tail = TRUE)
    ) |>
    dplyr::select(.data$subject, .data$n_obs, .data$Di, .data$Si, .data$mean_resid, .data$PIT,
                  .data$chisq_threshold, .data$significant, .data$interaction) |>
    dplyr::arrange(.data$subject)

  pit <- pmin(pmax(result_tbl$PIT, 1e-12), 1 - 1e-12)
  if (length(pit) >= 3 && length(unique(pit)) >= 3) {
    ad_res  <- goftest::ad.test(pit, null = "punif")
    ad_stat <- as.numeric(ad_res$statistic)
    ad_p    <- as.numeric(ad_res$p.value)
  } else { ad_stat <- NA_real_; ad_p <- NA_real_ }

  overall_interaction <- if (p_value_total >= alpha) "additive" else {
    sig_counts <- dplyr::filter(result_tbl, .data$significant) |>
      dplyr::count(.data$interaction) |>
      dplyr::filter(.data$interaction %in% c("synergistic","antagonistic"))
    if (nrow(sig_counts) == 0) "undetermined" else sig_counts$interaction[which.max(sig_counts$n)]
  }

  chi_square_result <- tibble::tibble(
    chi2_statistic = D_total,
    df             = df_total,
    p_value        = p_value_total,
    overall_interaction = overall_interaction,
    AD_statistic_goodness_of_fit = ad_stat,
    AD_p_value_goodness_of_fit   = ad_p
  )

  x_max <- max(D_total, stats::qchisq(0.999, df = df_total))
  x_seq <- seq(0, x_max, length.out = 500)
  curve_df <- tibble::tibble(x = x_seq, y = stats::dchisq(x_seq, df = df_total))
  shade_df <- dplyr::filter(curve_df, .data$x >= D_total)

  p <- ggplot2::ggplot(curve_df, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_line() +
    ggplot2::geom_area(data = shade_df, ggplot2::aes(x = .data$x, y = .data$y), alpha = 0.6) +
    ggplot2::geom_vline(xintercept = D_total, linetype = "dashed", linewidth = 1.2) +
    ggplot2::labs(
      title = "Chi-squared test for group-level interaction",
      subtitle = paste0("Statistic = ", round(D_total, 2),
                        ", df = ", df_total,
                        ", p = ", format(p_value_total, scientific = TRUE, digits = 3)),
      x = "Chi-squared statistic", y = "Density"
    ) + ggplot2::theme_bw()

  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    utils::write.csv(chi_square_result, file.path(out_dir, "chi_square_result.csv"), row.names = FALSE)
    utils::write.csv(result_tbl,        file.path(out_dir, "chi_square_subject_results.csv"), row.names = FALSE)
    ggplot2::ggsave(file.path(out_dir, "chi_square_test_plot.png"), p, width = 8, height = 6, dpi = 300)
  }

  list(overall = chi_square_result, per_subject = result_tbl, plot = p)
}
