# plot_individual.R
#' Plot subject-level observations, predictions, and additive curve
#' @param fit synpdx_fit
#' @param data data.frame with columns subject,time,treatment,logvol_norm
#' @param control control arm name
#' @param drug_a drug A arm name
#' @param drug_b drug B arm name
#' @param combo combination arm name
#' @param cutoff "auto" or "none" (auto = last common time)
#' @param out_dir optional directory to save PNG
#' @param width plot width
#' @param height plot height
#' @param dpi dpi
#' @export
synpdx_plot_individual <- function(fit, data, control, drug_a, drug_b, combo,
                                   cutoff = c("auto","none"), out_dir = NULL,
                                   width = 15, height = 12, dpi = 300) {
  cutoff <- match.arg(cutoff)
  stopifnot(inherits(fit, "synpdx_fit"))
  selected_model <- fit$selected_model
  data_processed <- fit$data_processed
  tmax_tbl <- dplyr::rename(fit$tmax_tbl, cutoff = .data$tmax)

  subjects_in_model_order <- unique(selected_model@data@data$subject)
  indiv_params <- as.data.frame(saemix::psi(selected_model))
  indiv_params$subject <- subjects_in_model_order
  indiv_params <- dplyr::relocate(indiv_params, .data$subject)

  observed_control_mono <- dplyr::filter(data_processed, .data$treatment != combo)
  observed_combination  <- dplyr::filter(data_processed, .data$treatment == combo)

  prediction_grid_base <- data_processed |>
    dplyr::group_by(.data$subject) |>
    dplyr::reframe(time = seq(0, max(.data$time), by = 0.5), .groups = "drop")
  grid_control <- dplyr::mutate(prediction_grid_base, I_A = 0L, I_B = 0L)
  grid_drug_a  <- dplyr::mutate(prediction_grid_base, I_A = 1L, I_B = 0L)
  grid_drug_b  <- dplyr::mutate(prediction_grid_base, I_A = 0L, I_B = 1L)

  model_func <- selected_model@model@model
  predict_on_grid <- function(grid, params) {
    data_to_predict <- dplyr::inner_join(grid, params, by = "subject")
    all_psi <- as.matrix(tibble::column_to_rownames(params, "subject"))
    all_subjects <- rownames(all_psi)
    dplyr::group_by(data_to_predict, .data$subject) |>
      dplyr::group_modify(function(df, key) {
        sid  <- which(all_subjects == key$subject[[1]])
        xidep_df <- df[, c("time","I_A","I_B")]
        df$predicted_logvol <- model_func(all_psi, rep(sid, nrow(df)), xidep_df)
        df
      }) |>
      dplyr::ungroup()
  }

  pred_control <- predict_on_grid(grid_control, indiv_params)
  pred_drug_a  <- predict_on_grid(grid_drug_a,  indiv_params)
  pred_drug_b  <- predict_on_grid(grid_drug_b,  indiv_params)

  pred_additive <- pred_drug_a |>
    dplyr::select(.data$subject, .data$time, pred_a = .data$predicted_logvol) |>
    dplyr::left_join(dplyr::select(pred_drug_b, .data$subject, .data$time, pred_b = .data$predicted_logvol),
                     by = c("subject","time")) |>
    dplyr::left_join(dplyr::select(pred_control, .data$subject, .data$time, pred_c = .data$predicted_logvol),
                     by = c("subject","time")) |>
    dplyr::mutate(predicted_logvol = .data$pred_a + .data$pred_b - .data$pred_c) |>
    dplyr::select(.data$subject, .data$time, .data$predicted_logvol)

  tidy_predictions <- dplyr::bind_rows(
    dplyr::mutate(pred_control, curve = control),
    dplyr::mutate(pred_drug_a,  curve = drug_a),
    dplyr::mutate(pred_drug_b,  curve = drug_b),
    dplyr::mutate(pred_additive, curve = "additive")
  ) |>
    dplyr::left_join(tmax_tbl, by = "subject") |>
    dplyr::mutate(line_type = ifelse(.data$time <= .data$cutoff, "solid", "dashed"))

  color_map <- c(control = "#1b9e77", additive = "black", combination = "black", prediction_cutoff = "red")
  names(color_map)[1] <- control
  color_map[drug_a] <- "#d95f02"
  color_map[drug_b] <- "#7570b3"
  breaks <- c(control, drug_a, drug_b, combo, "additive")
  labels <- c(paste(control,"(Obs.)"), paste(drug_a,"(Obs.)"), paste(drug_b,"(Obs.)"),
              paste(combo,"(Obs.)"), "Additive (Pred.)")

  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = observed_control_mono,
                        ggplot2::aes(x = .data$time, y = .data$logvol_norm, color = .data$treatment), size = 2) +
    ggplot2::geom_point(data = observed_combination,
                        ggplot2::aes(x = .data$time, y = .data$logvol_norm, color = combo),
                        shape = 1, size = 2, stroke = 1) +
    ggplot2::geom_line(data = dplyr::filter(tidy_predictions, .data$curve != "additive", .data$line_type == "solid"),
                       ggplot2::aes(x = .data$time, y = .data$predicted_logvol, color = .data$curve,
                                    group = interaction(.data$subject, .data$curve)), linetype = "solid", linewidth = 0.8) +
    ggplot2::geom_line(data = dplyr::filter(tidy_predictions, .data$curve != "additive", .data$line_type == "dashed"),
                       ggplot2::aes(x = .data$time, y = .data$predicted_logvol, color = .data$curve,
                                    group = interaction(.data$subject, .data$curve)), linetype = "dashed", linewidth = 0.8) +
    ggplot2::geom_line(data = dplyr::filter(tidy_predictions, .data$curve == "additive", .data$line_type == "solid"),
                       ggplot2::aes(x = .data$time, y = .data$predicted_logvol, color = "additive", group = .data$subject),
                       linetype = "solid", linewidth = 0.9) +
    ggplot2::geom_line(data = dplyr::filter(tidy_predictions, .data$curve == "additive", .data$line_type == "dashed"),
                       ggplot2::aes(x = .data$time, y = .data$predicted_logvol, color = "additive", group = .data$subject),
                       linetype = "dashed", linewidth = 0.9) +
    ggplot2::geom_vline(data = tmax_tbl, ggplot2::aes(xintercept = .data$cutoff, linetype = "prediction_cutoff"),
                        color = "red", linewidth = 0.8) +
    ggh4x::facet_wrap2(~ subject, scales = "free_y", axes = "all") +
    ggplot2::scale_color_manual(name = NULL, values = color_map, breaks = breaks, labels = labels) +
    ggplot2::scale_linetype_manual(name = NULL,
                                   values = c("solid"="solid","dashed"="dashed","prediction_cutoff"="dotted"),
                                   labels = c("solid"="Predicted","dashed"="Extrapolated","prediction_cutoff"="Prediction Cutoff")) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(shape = c(19,19,19,1,NA),
                                                                      linetype = c("blank","blank","blank","blank","solid")),
                                                  nrow = 2, order = 1),
                    linetype = ggplot2::guide_legend(nrow = 1, order = 2)) +
    ggplot2::labs(x = "Time (days)", y = "Normalized Log-Volume") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "bottom",
                   strip.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(colour = "black", linewidth = 1.0))

  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(out_dir, "individual_prediction_plot.png"),
                    plot = p, width = width, height = height, dpi = dpi)
  }
  p
}
