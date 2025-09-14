#' Fit and analyze synpdx in one call
#'
#' Runs model selection + fitting, plots subject predictions, computes per-subject
#' distances, performs chi-square test and bootstrap CIs calculation for hypothesis testing.
#'
#' @param data data.frame with columns \code{subject}, \code{time}, \code{treatment}, \code{logvol_norm}.
#' @param control character. Control arm name.
#' @param drug_a character. Drug A arm name.
#' @param drug_b character. Drug B arm name.
#' @param combo  character. Combination arm name.
#' @param B integer. Bootstrap replicates.
#' @param out_dir character or NULL. Output directory for figures/CSVs.
#' @param seed integer. Random seed for SAEM fitting and bootstrap resampling.
#'
#' @return list(fit, plot, Di, chisq, bootstrap)
#'
#' @examples
#' \dontrun{
#' # Example 1: run on your CSV (Windows path -> use forward slashes)
#' library(synpdx)
#' csv <- "C:/Users/data.csv"
#' df  <- read.csv(csv, stringsAsFactors = FALSE)
#'
#' # Required columns: subject, time, treatment, logvol_norm
#' stopifnot(all(c("subject","time","treatment","logvol_norm") %in% names(df)))
#'
#' res <- synpdx_fit(
#'   df,
#'   control = "control",
#'   drug_a  = "BYL719",
#'   drug_b  = "binimetinib",
#'   combo   = "combination",
#'   B       = 200,
#'   seed    = 1234,
#'   out_dir = "synpdx_outputs"
#' )
#' print(res$chisq$overall)
#' res$bootstrap$summary
#' res$plot
#'
#' # Example 2: if your file has raw volumes 'volume' instead of log-normalized:
#' # Create logvol_norm = log(V(t)/V(0)) per subject-arm.
#' # library(dplyr)
#' # df <- df %>%
#' #   dplyr::group_by(subject, treatment) %>%
#' #   dplyr::arrange(time, .by_group = TRUE) %>%
#' #   dplyr::mutate(logvol_norm = log(volume / dplyr::first(volume))) %>%
#' #   dplyr::ungroup()
#' }
#'
#' @export
synpdx_fit <- function(data, control, drug_a, drug_b, combo,
                       B = 100, out_dir = NULL, seed = 990202) {
  set.seed(seed)
  fit <- synpdx_model_fit(data, control, drug_a, drug_b, combo, seed = seed)
  p   <- synpdx_plot_individual(fit, data, control, drug_a, drug_b, combo, out_dir = out_dir)
  Di  <- synpdx_compute_Di(fit, data, control, drug_a, drug_b, combo, out_dir = out_dir)
  cs  <- synpdx_chisq_test(Di$Di_vec, Di$combo_df, out_dir = out_dir)
  bt  <- synpdx_bootstrap(fit, data, control, drug_a, drug_b, combo,
                          B = B, out_dir = out_dir, seed = seed)
  list(fit = fit, plot = p, Di = Di, chisq = cs, bootstrap = bt)
}
