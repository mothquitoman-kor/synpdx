# R/fit.R

# ---- robust BIC extractor for SaemixObject -----------------------------------
get_bic_saemix <- function(obj) {
  val <- try(saemix::AIC(obj), silent = TRUE)
  if (!inherits(val, "try-error")) {
    nms <- names(val)
    if (!is.null(nms)) {
      cand_names <- c("BIC.is","BIC.IS","BIC","BIC.lin")
      for (nm in cand_names) if (nm %in% nms && is.finite(val[[nm]])) return(unname(val[[nm]]))
      idx <- grep("BIC", nms, ignore.case = TRUE)
      if (length(idx)) {
        v <- as.numeric(val[idx[1]])
        if (is.finite(v)) return(v)
      }
    } else {
      v <- as.numeric(val)[1]
      if (is.finite(v)) return(v)
    }
  }
  ll <- try(obj@results@ll.is,  silent = TRUE)
  if (!(is.numeric(ll) && length(ll)==1 && is.finite(ll))) ll <- try(obj@results@ll.lin, silent = TRUE)
  if (is.numeric(ll) && length(ll)==1 && is.finite(ll)) {
    k <- try(length(obj@results@fixed.effects), silent = TRUE); if (!is.numeric(k) || !is.finite(k)) k <- 0
    n <- try(nrow(obj@data@data),              silent = TRUE); if (!is.numeric(n) || !is.finite(n)) n <- 0
    if (n > 0) return(-2 * ll + k * log(n))
  }
  NA_real_
}

# ---- main fitter --------------------------------------------------------------
#' Internal: Fit one or auto-select NLME model(s) with optional RE backward selection
#' @keywords internal
synpdx_model_fit <- function(dat, control, drug_a, drug_b, combo,
                             seed = 990202,
                             model = c("auto","exp_const","exp_decay","gomp_const","gomp_decay","logistic_const","logistic_decay"),
                             select_random_effects = FALSE) {  # default FALSE to avoid unintended selection
  model <- match.arg(model, c("auto","exp_const","exp_decay","gomp_const","gomp_decay","logistic_const","logistic_decay"))
  set.seed(seed)

  # ---------- checks ----------
  dat <- as.data.frame(dat)
  if (!nrow(dat)) stop("`dat` has 0 rows")
  need <- c("subject","time","treatment","logvol_norm")
  miss <- setdiff(need, names(dat))
  if (length(miss)) stop("`dat` missing: ", paste(miss, collapse = ", "))

  # ---------- indicators, keep all 4 arms ----------
  data_processed <- dat
  data_processed$subject <- as.character(data_processed$subject)
  data_processed$I_A <- as.integer(data_processed$treatment %in% c(drug_a, combo))
  data_processed$I_B <- as.integer(data_processed$treatment %in% c(drug_b, combo))
  keep_arms <- c(control, drug_a, drug_b, combo)
  data_processed <- data_processed[data_processed$treatment %in% keep_arms, , drop = FALSE]
  data_processed <- data_processed[order(data_processed$subject, data_processed$time), , drop = FALSE]

  # subject-wise common last time across 4 arms (for plotting/cutoff)
  tmax_list <- lapply(split(data_processed, data_processed$subject), function(subd){
    gmax <- tapply(subd$time, subd$treatment, max, na.rm = TRUE)
    tcut <- if (all(keep_arms %in% names(gmax))) min(gmax[keep_arms]) else NA_real_
    data.frame(subject = as.character(subd$subject[1]), tmax = tcut, stringsAsFactors = FALSE)
  })
  tmax_tbl <- if (length(tmax_list)) do.call(rbind, tmax_list) else data.frame(subject = character(), tmax = numeric())

  # ---------- training data (control + A + B) ----------
  df <- subset(data_processed, treatment %in% c(control, drug_a, drug_b))
  if (!nrow(df)) stop("No rows for control/monotherapies")
  df <- df[order(df$subject, df$time), , drop = FALSE]

  saem_data <- saemix::saemixData(
    name.data       = df,
    name.group      = "subject",
    name.predictors = c("time","I_A","I_B"),
    name.response   = "logvol_norm"
  )

  # ---------- model library (defined in models.R) ----------
  model_map <- list(
    exp_const      = exp_const,
    exp_decay      = exp_decay,
    gomp_const     = gomp_const,
    gomp_decay     = gomp_decay,
    logistic_const = logistic_const,
    logistic_decay = logistic_decay
  )

  # ---------- fit single model with/without RE selection ----------
  fit_one <- function(modfun, do_re_selection) {
    th0    <- attr(modfun, "theta0")
    pnames <- attr(modfun, "parnames")
    if (is.null(th0) || is.null(pnames)) stop("theta0/parnames missing in model")
    if (length(th0) != length(pnames)) stop("length(theta0) != length(parnames)")

    psi0 <- matrix(th0, ncol = length(th0), byrow = TRUE)
    colnames(psi0) <- pnames

    # default: diagonal random effects (no selection)
    cov_use <- diag(1, length(pnames))
    dimnames(cov_use) <- list(pnames, pnames)

    if (isTRUE(do_re_selection)) {
      if (!exists("backward_selection_saemix", mode = "function"))
        stop("backward_selection_saemix() not found (required when select_random_effects=TRUE)")
      sel <- backward_selection_saemix(
        par_names        = pnames,
        saemix_data      = saem_data,
        saemix_model_fun = modfun,
        psi0             = psi0,
        alpha            = 0.01
      )
      if (!is.null(sel$covariance.model)) cov_use <- sel$covariance.model
    }

    saem_model <- saemix::saemixModel(
      model            = modfun,
      psi0             = psi0,
      transform.par    = rep(0, length(th0)),
      error.model      = "constant",
      name.modpar      = pnames,
      covariance.model = cov_use
    )
    ctrl <- saemix::saemixControl(seed = seed, displayProgress = FALSE)
    saemix::saemix(saem_model, saem_data, control = ctrl)
  }

  # ---------- auto selection over 6 candidates by BIC ----------
  if (model == "auto") {
    cand <- names(model_map)
    fits <- vector("list", length(cand))
    bics <- rep(Inf, length(cand))

    for (i in seq_along(cand)) {
      mf  <- model_map[[cand[i]]]
      obj <- try(fit_one(mf, do_re_selection = isTRUE(select_random_effects)), silent = TRUE)
      if (inherits(obj, "SaemixObject")) {
        fits[[i]] <- obj
        bic_i <- suppressWarnings(get_bic_saemix(obj))
        if (is.finite(bic_i)) bics[i] <- bic_i
      }
    }
    if (!any(is.finite(bics)))
      stop("auto selection failed: BIC not available from any fitted model")

    pick_idx <- which.min(bics)
    saem_obj <- fits[[pick_idx]]
    picked   <- cand[pick_idx]

  } else {
    if (!model %in% names(model_map)) stop("Unknown model: ", model)
    saem_obj <- fit_one(model_map[[model]], do_re_selection = isTRUE(select_random_effects))
    picked   <- model
  }

  # ---------- return ----------
  out <- list(
    saemix               = saem_obj,
    selected_model       = saem_obj,
    selected_model_name  = picked,
    data_processed       = data_processed,
    tmax_tbl             = tmax_tbl
  )
  class(out) <- c("synpdx_fit", class(out))
  out
}
