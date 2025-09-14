#' @importFrom rlang .data
NULL

# silence NSE notes for dplyr/ggplot2 columns
utils::globalVariables(c(
  ".", ".data", "subject", "treatment", "time", "tmax", "tmax_arm", "arm",
  "logvol_norm", "predicted_logvol", "pred_a", "pred_b", "pred_c",
  "line_type", "curve", "DV", "Di", "Si", "PIT", "n_obs",
  "chisq_threshold", "mean_resid", "sign_resid", "significant",
  "x", "y", "density", "cutoff", "I_A", "I_B", "after_stat"
))
