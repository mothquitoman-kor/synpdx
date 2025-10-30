# backward_selection.R
#' Boundary LRT backward selection for random effects (saemix)
#' @param par_names character vector of parameter names (log-scale)
#' @param saemix_data saemixData
#' @param saemix_model_fun function(psi,id,xidep) returning y
#' @param psi0 1xP matrix of initial fixed effects (log-scale)
#' @param alpha boundary LRT alpha, default 0.01
#' @return list(selected, covariance.model)
#' @keywords internal
backward_selection_saemix <- function(par_names, saemix_data, saemix_model_fun, psi0, alpha = 0.01) {
  stopifnot(length(par_names) == ncol(psi0))
  P <- length(par_names)
  keep <- par_names
  mask_from <- function(keep) {
    M <- matrix(0, P, P, dimnames = list(par_names, par_names))
    if (length(keep)) M[keep, keep] <- diag(1, length(keep))
    M
  }
  full_model <- saemix::saemixModel(
    model = saemix_model_fun, psi0 = psi0,
    transform.par = rep(0, P), covariance.model = mask_from(keep),
    error.model = "constant"
  )
  ctrl <- saemix::saemixControl(fim = FALSE, map = TRUE, ll.is = TRUE,
                                nbiter.saemix = c(300, 100), nb.chains = 2,
                                print = FALSE, displayProgress = FALSE, save.graphs = FALSE)
  fit_full <- saemix::saemix(full_model, saemix_data, control = ctrl)
  ll_full  <- as.numeric(fit_full@results@ll.is)

  improved <- TRUE
  while (improved && length(keep) > 0) {
    improved <- FALSE
    candidates <- keep
    pvals <- rep(NA_real_, length(candidates)); names(pvals) <- candidates
    for (rmv in candidates) {
      keep_try <- setdiff(keep, rmv)
      model_try <- saemix::saemixModel(
        model = saemix_model_fun, psi0 = psi0,
        transform.par = rep(0, P), covariance.model = mask_from(keep_try),
        error.model = "constant"
      )
      fit_try <- try(saemix::saemix(model_try, saemix_data, control = ctrl), silent = TRUE)
      if (inherits(fit_try, "try-error")) next
      ll_try <- as.numeric(fit_try@results@ll.is)
      LR <- 2 * (ll_full - ll_try)
      pvals[rmv] <- 0.5 * (1 - stats::pchisq(LR, df = 1))
    }
    if (all(is.na(pvals))) break
    worst <- names(which.max(pvals))
    if (!is.na(pvals[worst]) && pvals[worst] > alpha) {
      keep <- setdiff(keep, worst)
      full_model <- saemix::saemixModel(
        model = saemix_model_fun, psi0 = psi0,
        transform.par = rep(0, P), covariance.model = mask_from(keep),
        error.model = "constant"
      )
      fit_full <- saemix::saemix(full_model, saemix_data, control = ctrl)
      ll_full  <- as.numeric(fit_full@results@ll.is)
      improved <- TRUE
    }
  }
  list(selected = keep, covariance.model = mask_from(keep))
}
