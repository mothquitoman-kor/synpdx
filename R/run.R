# R/run.R
#' synpdx: NLME-based interaction analysis for PDX combinations
#' @keywords internal
"_PACKAGE"

#' Fit NLME models and run interaction tests
#'
#' This high-level wrapper performs:
#' 1) fitting of an NLME tumor-growth and drug-effect model on control and monotherapy data,
#' 2) prediction and distance computation (D_i) for the combination arm,
#' 3) chi-square additivity testing, and
#' 4) bootstrap testing of standardized signed distances.
#'
#' Input data must already include a log-scale response column named `logvol_norm`.
#'
#' @param dat A data.frame containing:
#'   - subject: character or coercible identifier for each PDX model
#'   - time: numeric observation time
#'   - treatment: character treatment arm label (must contain the four arm names below)
#'   - logvol_norm: numeric log-transformed tumor volume
#'   Indicators I_A and I_B are computed internally.
#'
#' @param control Character. Label of the control arm in `treatment` (e.g., "control").
#' @param drug_a Character. Label of the monotherapy A arm (e.g., "BYL719").
#' @param drug_b Character. Label of the monotherapy B arm (e.g., "binimetinib").
#' @param combo Character. Label of the combination arm (e.g., "combination").
#' @param B Integer ≥1. Number of bootstrap replicates. Default is 100.
#' @param out_dir NULL or character path. If not NULL, CSV and PNG outputs are written there.
#' @param seed Integer random seed. Default is 990202.
#' @param cutoff One of c("auto","none").
#'   - "auto": predictions in the combination arm are truncated at the smallest last-time point
#'     shared across all four arms (control, A, B, combination).
#'   - "none": use all available time points.
#' @param model One of c("auto","exp_const","exp_decay","gomp_const","gomp_decay","logistic_const","logistic_decay").
#'   - "auto": all six models are fitted and compared by BIC; the one with lowest BIC is chosen.
#'   - Otherwise, the specified model is used directly.
#' @param select_random_effects Logical.
#'   - TRUE: perform backward likelihood ratio tests for random-effect selection using backward_selection_saemix.
#'   - FALSE: skip selection and use a diagonal covariance (one random effect per parameter).
#'   In either case, if model="auto", the final decision among six models is based on BIC.
#'
#' @return A list with components:
#'   - fit: synpdx_fit object (selected model and saemix result)
#'   - plot: ggplot object from synpdx_plot_individual
#'   - Di: list from synpdx_compute_Di
#'   - chisq: list from synpdx_chisq_test
#'   - bootstrap: list from synpdx_bootstrap
#'
#' @examples
#' \dontrun{
#' # Example dataset provided with the package
#' csv_path <- system.file(
#'   "extdata", "PDX_CRC_BYL719+binimetinib_curated.csv",
#'   package = "synpdx"
#' )
#' df <- read.csv(csv_path, stringsAsFactors = FALSE)
#'
#' # ex1: auto model selection, cutoff="auto", random effects selection=TRUE
#' res1 <- synpdx_fit(
#'   dat = df,
#'   control = "control", drug_a = "BYL719", drug_b = "binimetinib", combo = "combination",
#'   B = 100, seed = 1, cutoff = "auto",
#'   model = "auto",
#'   select_random_effects = TRUE,
#'   out_dir = "synpdx_out_ex1"
#' )
#'
#' # ex2: gomp_decay model, cutoff="none", random effects selection=FALSE
#' res2 <- synpdx_fit(
#'   dat = df,
#'   control = "control", drug_a = "BYL719", drug_b = "binimetinib", combo = "combination",
#'   B = 100, seed = 2, cutoff = "none",
#'   model = "gomp_decay",
#'   select_random_effects = FALSE,
#'   out_dir = "synpdx_out_ex2"
#' )
#'
#' # ex3: logistic_const model, cutoff="auto", random effects selection=TRUE
#' res3 <- synpdx_fit(
#'   dat = df,
#'   control = "control", drug_a = "BYL719", drug_b = "binimetinib", combo = "combination",
#'   B = 50, seed = 3, cutoff = "auto",
#'   model = "logistic_const",
#'   select_random_effects = TRUE,
#'   out_dir = "synpdx_out_ex3"
#' )
#'
#' # ex4: exp_decay model, cutoff="none", random effects selection=FALSE
#' res4 <- synpdx_fit(
#'   dat = df,
#'   control = "control", drug_a = "BYL719", drug_b = "binimetinib", combo = "combination",
#'   B = 30, seed = 4, cutoff = "none",
#'   model = "exp_decay",
#'   select_random_effects = FALSE,
#'   out_dir = "synpdx_out_ex4"
#' )
#' }
#'
#' @export
synpdx_fit <- function(dat, control, drug_a, drug_b, combo,
                       B = 100, out_dir = NULL, seed = 990202,
                       cutoff = c("auto","none"),
                       model = c("auto","exp_const","exp_decay","gomp_const","gomp_decay","logistic_const","logistic_decay"),
                       select_random_effects = TRUE) {

  cutoff <- match.arg(cutoff)
  model  <- match.arg(model)

  if (is.function(dat)) {
    nm <- deparse(substitute(dat))
    stop("`dat` resolved to a function (", nm, "). Pass a data.frame via `dat=`.")
  }
  dat <- as.data.frame(dat)
  if (!nrow(dat)) stop("`dat` has 0 rows")

  fit <- synpdx_model_fit(
    dat = dat, control = control, drug_a = drug_a, drug_b = drug_b, combo = combo,
    seed = seed, model = model, select_random_effects = select_random_effects
  )

  p  <- synpdx_plot_individual(
    fit, dat, control, drug_a, drug_b, combo,
    cutoff = cutoff, out_dir = out_dir
  )

  Di <- synpdx_compute_Di(
    fit, dat, control, drug_a, drug_b, combo,
    cutoff = cutoff, out_dir = out_dir
  )

  cs <- synpdx_chisq_test(
    Di_vec = Di$Di_vec, combo_df = Di$combo_df, out_dir = out_dir
  )

  bt <- synpdx_bootstrap(
    fit, dat, control, drug_a, drug_b, combo,
    B = B, seed = seed, cutoff = cutoff, out_dir = out_dir
  )

  list(fit = fit, plot = p, Di = Di, chisq = cs, bootstrap = bt)
}
