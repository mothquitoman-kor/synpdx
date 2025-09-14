#' Internal: fit NLME growth–drug model(s) with SAEM and select by BIC
#' @keywords internal
synpdx_model_fit <- function(data, control, drug_a, drug_b, combo,
                             alpha_eta = 0.01, seed = 990202) {
  stopifnot(all(c("subject","time","treatment","logvol_norm") %in% names(data)))
  data <- dplyr::mutate(data, subject = as.character(.data$subject))
  arms_all <- c(control, drug_a, drug_b, combo)

  data_processed <- data |>
    dplyr::mutate(
      I_A = as.integer(.data$treatment == drug_a | .data$treatment == combo),
      I_B = as.integer(.data$treatment == drug_b | .data$treatment == combo)
    )

  data_for_fitting <- dplyr::filter(data_processed, .data$treatment %in% c(control, drug_a, drug_b))

  sdm <- saemix::saemixData(
    name.data = data_for_fitting,
    name.group = "subject",
    name.predictors = c("time","I_A","I_B"),
    name.response   = "logvol_norm",
    name.X = "time"
  )

  tmax_tbl <- data_processed |>
    dplyr::filter(.data$treatment %in% arms_all) |>
    dplyr::group_by(.data$subject, .data$treatment) |>
    dplyr::summarise(tmax_arm = max(.data$time), .groups = "drop") |>
    dplyr::group_by(.data$subject) |>
    dplyr::summarise(tmax = min(.data$tmax_arm), .groups = "drop")

  model_list <- list(
    exp_const       = list(fun = exp_const,       eta = c("l_a","l_beta_a","l_beta_b")),
    exp_decay       = list(fun = exp_decay,       eta = c("l_a","l_beta_a","l_beta_b","l_k_a","l_k_b")),
    gomp_const      = list(fun = gomp_const,      eta = c("l_alpha1","l_alpha2","l_beta_a","l_beta_b")),
    gomp_decay      = list(fun = gomp_decay,      eta = c("l_alpha1","l_alpha2","l_beta_a","l_beta_b","l_k_a","l_k_b")),
    logistic_const  = list(fun = logistic_const,  eta = c("l_K","l_r","l_beta_a","l_beta_b")),
    logistic_decay  = list(fun = logistic_decay,  eta = c("l_K","l_r","l_beta_a","l_beta_b","l_k_a","l_k_b"))
  )

  psi0_list <- list(
    exp_const       = matrix(log(c(a=0.2, beta_a=0.01, beta_b=0.01)),
                             nrow=1, dimnames=list(NULL, c("l_a","l_beta_a","l_beta_b"))),
    exp_decay       = matrix(log(c(a=0.2, beta_a=0.01, beta_b=0.01, k_a=0.05, k_b=0.05)),
                             nrow=1, dimnames=list(NULL, c("l_a","l_beta_a","l_beta_b","l_k_a","l_k_b"))),
    gomp_const      = matrix(log(c(alpha1=0.2, alpha2=0.01, beta_a=0.01, beta_b=0.01)),
                             nrow=1, dimnames=list(NULL, c("l_alpha1","l_alpha2","l_beta_a","l_beta_b"))),
    gomp_decay      = matrix(log(c(alpha1=0.2, alpha2=0.01, beta_a=0.01, beta_b=0.01, k_a=0.05, k_b=0.05)),
                             nrow=1, dimnames=list(NULL, c("l_alpha1","l_alpha2","l_beta_a","l_beta_b","l_k_a","l_k_b"))),
    logistic_const  = matrix(log(c(K=4, r=0.1, beta_a=0.01, beta_b=0.01)),
                             nrow=1, dimnames=list(NULL, c("l_K","l_r","l_beta_a","l_beta_b"))),
    logistic_decay  = matrix(log(c(K=4, r=0.1, beta_a=0.01, beta_b=0.01, k_a=0.05, k_b=0.05)),
                             nrow=1, dimnames=list(NULL, c("l_K","l_r","l_beta_a","l_beta_b","l_k_a","l_k_b")))
  )

  set.seed(seed)
  ctrl <- saemix::saemixControl(fim = TRUE, map = TRUE, ll.is = TRUE,
                                nbiter.saemix = c(400, 150), nb.chains = 3,
                                fix.seed = TRUE, seed = seed,
                                displayProgress = FALSE, print = FALSE, save.graphs = FALSE)

  final_fits <- list()
  for (nm in names(model_list)) {
    m    <- model_list[[nm]]
    psi0 <- psi0_list[[nm]]
    bw <- backward_selection_saemix(
      par_names = m$eta, saemix_data = sdm, saemix_model_fun = m$fun, psi0 = psi0, alpha = alpha_eta
    )
    model_final <- saemix::saemixModel(
      model = m$fun, psi0 = psi0, transform.par = rep(0, length(m$eta)),
      covariance.model = bw$covariance.model, error.model = "constant"
    )
    final_fits[[nm]] <- saemix::saemix(model_final, sdm, control = ctrl)
  }

  bic_values <- vapply(final_fits, function(f) {
    val <- try(stats::BIC(f, method = "is"), silent = TRUE)
    if (inherits(val,"try-error") || is.na(val)) stats::BIC(f, method = "lin") else val
  }, numeric(1))

  best_model_name <- names(bic_values)[which.min(bic_values)]
  selected_model  <- final_fits[[best_model_name]]

  obj <- list(
    selected_model = selected_model,
    best_model_name = best_model_name,
    tmax_tbl = tmax_tbl,
    data_processed = data_processed,
    arms = list(control=control, drug_a=drug_a, drug_b=drug_b, combo=combo)
  )
  class(obj) <- "synpdx_fit"
  obj
}

#' Print a synpdx_fit summary
#' @param x synpdx_fit
#' @param ... unused
#' @export
print.synpdx_fit <- function(x, ...) {
  cat("synpdx_fit: best model =", x$best_model_name, "\n")
  invisible(x)
}
