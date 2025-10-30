# di.R
#' Compute per-subject distance D_i and combo predictions
#' @export
synpdx_compute_Di <- function(fit, data, control, drug_a, drug_b, combo,
                              cutoff = c("auto", "none"), out_dir = NULL) {
  cutoff <- match.arg(cutoff)
  grDevices::devAskNewPage(FALSE)

  # required columns
  stopifnot(all(c("subject","time","treatment","logvol_norm") %in% names(data)))

  # base df with indicators
  df <- data
  df <- dplyr::mutate(
    df,
    subject = as.character(.data$subject),
    I_A = as.integer(.data$treatment %in% c(drug_a, combo)),
    I_B = as.integer(.data$treatment %in% c(drug_b, combo))
  )

  # locate SaemixObject
  get_saem <- function(x) {
    if (inherits(x, "SaemixObject")) return(x)
    cand <- list(
      tryCatch(x$saemix, error = function(...) NULL),
      tryCatch(x$fit,    error = function(...) NULL),
      tryCatch(x$object, error = function(...) NULL),
      tryCatch(x$saem,   error = function(...) NULL),
      tryCatch(if (!is.null(x$fit)) x$fit$saemix else NULL, error = function(...) NULL)
    )
    for (y in cand) if (inherits(y, "SaemixObject")) return(y)
    if (is.list(x)) {
      for (nm in names(x)) {
        y <- x[[nm]]
        if (inherits(y, "SaemixObject")) return(y)
        if (is.list(y)) {
          for (nm2 in names(y)) {
            z <- y[[nm2]]
            if (inherits(z, "SaemixObject")) return(z)
          }
        }
      }
    }
    stop("Could not locate SaemixObject inside `fit`.")
  }
  saem <- get_saem(fit)

  # model function
  model_fun <- NULL
  try(model_fun <- saem@model@model, silent = TRUE)
  if (!is.function(model_fun)) stop("Model function not found in SaemixObject.")

  # ---- parameter names from model as single source of truth ----
  pnames <- try(saem@model@name.modpar, silent = TRUE)
  if (inherits(pnames, "try-error") || is.null(pnames) || !length(pnames))
    stop("Parameter names not found in SaemixObject")

  # ---- fixed effects as named vector aligned to pnames ----
  mu_raw <- try(saem@results@fixed.effects, silent = TRUE)
  if (inherits(mu_raw, "try-error") || is.null(mu_raw) || length(mu_raw) == 0)
    stop("Fixed effects not found")

  if (is.matrix(mu_raw)) {
    rn <- rownames(mu_raw)
    if (!is.null(rn) && all(pnames %in% rn)) {
      mu_full <- stats::setNames(as.numeric(mu_raw[pnames, 1]), pnames)
    } else {
      stopifnot(nrow(mu_raw) >= length(pnames))
      mu_full <- stats::setNames(as.numeric(mu_raw[seq_along(pnames), 1]), pnames)
    }
  } else {
    nm <- names(mu_raw)
    if (!is.null(nm) && all(pnames %in% nm)) {
      mu_full <- stats::setNames(as.numeric(mu_raw[pnames]), pnames)
    } else {
      stopifnot(length(mu_raw) >= length(pnames))
      mu_full <- stats::setNames(as.numeric(mu_raw[seq_along(pnames)]), pnames)
    }
  }
  stopifnot(all(is.finite(mu_full)))

  # ---- EBEs from phi aligned to pnames ----
  get_ebes_df <- function(saem, df, pnames, mu_full) {
    phi <- try(as.data.frame(saem@results@phi), silent = TRUE)
    if (inherits(phi, "try-error") || is.null(phi) || ncol(phi) == 0)
      stop("Individual parameters (phi) unavailable")

    subcol <- if ("subject" %in% names(phi)) "subject" else {
      idcand <- names(phi)[match("id", tolower(names(phi)))]
      if (!is.na(idcand)) idcand else NULL
    }
    if (is.null(subcol)) { phi$subject <- sort(unique(df$subject))[seq_len(nrow(phi))]; subcol <- "subject" }
    phi$subject <- as.character(phi[[subcol]])

    cols <- match(pnames, names(phi))
    if (any(is.na(cols))) {
      num_cols <- which(vapply(phi, is.numeric, TRUE))
      if (subcol %in% names(phi)) num_cols <- setdiff(num_cols, match(subcol, seq_along(phi)))
      if (length(num_cols) < length(pnames)) stop("phi columns < parameter count")
      cols <- num_cols[seq_along(pnames)]
    }
    ebes_mat <- sweep(as.matrix(phi[, cols, drop = FALSE]), 2, mu_full[pnames], FUN = "-")
    ebes <- cbind(subject = phi$subject, as.data.frame(ebes_mat, check.names = FALSE))
    names(ebes) <- c("subject", pnames)
    ebes$subject <- as.character(ebes$subject)
    ebes
  }
  ebes <- get_ebes_df(saem, df, pnames, mu_full)

  # vectorized accessor
  mu_vec <- mu_full
  get_psi_for_subject <- function(sid) {
    idx <- match(sid, ebes$subject)
    if (is.na(idx)) stop("Subject ", sid, " not in EBEs")
    eta_i <- stats::setNames(as.numeric(ebes[idx, pnames, drop = TRUE]), pnames)
    mu_vec + eta_i
  }

  # ---- residual variance from training arms ----
  df_train <- dplyr::filter(df, .data$treatment %in% c(control, drug_a, drug_b))
  if (!nrow(df_train)) stop("No rows for control/monotherapies")
  train_res <- lapply(split(df_train, df_train$subject), function(subd) {
    psi_i <- get_psi_for_subject(unique(subd$subject))
    pred <- model_fun(
      psi   = matrix(psi_i, nrow = 1),
      id    = 1L,
      xidep = data.frame(time = subd$time, I_A = subd$I_A, I_B = subd$I_B)
    )
    as.numeric(pred) - subd$logvol_norm
  })
  resid_vec <- unlist(train_res, use.names = FALSE)
  sigma2_used <- mean(resid_vec^2, na.rm = TRUE)
  if (!is.finite(sigma2_used) || sigma2_used <= 0)
    sigma2_used <- stats::var(df_train$logvol_norm, na.rm = TRUE)

  # ---- combo predictions with cutoff handling ----
  df_combo <- dplyr::filter(df, .data$treatment == combo)
  if (!nrow(df_combo)) stop("No rows for combination arm")

  if (cutoff == "auto") {
    cut_tbl <- lapply(split(df, df$subject), function(subd) {
      gmax <- tapply(subd$time, subd$treatment, max, na.rm = TRUE)
      arms <- c(control, drug_a, drug_b, combo)
      if (!all(arms %in% names(gmax)))
        return(data.frame(subject = unique(subd$subject), t_cut = NA_real_))
      data.frame(subject = unique(subd$subject), t_cut = min(gmax[arms]))
    }) |> dplyr::bind_rows()
    df_combo <- dplyr::left_join(df_combo, cut_tbl, by = "subject") |>
      dplyr::filter(is.na(.data$t_cut) | .data$time <= .data$t_cut) |>
      dplyr::mutate(cutoff = .data$t_cut) |>
      dplyr::select(-.data$t_cut)
  } else {
    df_combo <- dplyr::mutate(df_combo, cutoff = NA_real_)
  }

  pred_list <- lapply(split(df_combo, df_combo$subject), function(subd) {
    psi_i <- get_psi_for_subject(unique(subd$subject))
    pred  <- model_fun(
      psi   = matrix(psi_i, nrow = 1),
      id    = 1L,
      xidep = data.frame(time = subd$time, I_A = 1L, I_B = 1L)
    )
    data.frame(
      subject = subd$subject,
      time = subd$time,
      DV = subd$logvol_norm,
      predicted_logvol = as.numeric(pred),
      I_A = 1L, I_B = 1L,
      cutoff = subd$cutoff[1],
      stringsAsFactors = FALSE
    )
  })
  combo_df <- dplyr::bind_rows(pred_list)

  # ---- D_i ----
  combo_df <- dplyr::mutate(combo_df, resid0 = .data$predicted_logvol - .data$DV)
  Di_tbl <- combo_df |>
    dplyr::group_by(.data$subject) |>
    dplyr::summarise(
      n_i = dplyr::n(),
      Di  = sum((.data$resid0^2) / sigma2_used, na.rm = TRUE),
      .groups = "drop"
    )
  Di_vec <- stats::setNames(Di_tbl$Di, Di_tbl$subject)

  # ---- outputs ----
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    utils::write.csv(combo_df[, c("subject","time","DV","predicted_logvol","I_A","I_B","cutoff")],
                     file.path(out_dir, "combo_pred_vs_obs.csv"), row.names = FALSE)
    utils::write.csv(data.frame(subject = names(Di_vec), Di = as.numeric(Di_vec)),
                     file.path(out_dir, "Di_by_subject.csv"), row.names = FALSE)
  }

  list(
    Di_vec = Di_vec,
    combo_df = combo_df[, c("subject","time","DV","predicted_logvol","I_A","I_B","cutoff")],
    sigma2_used = sigma2_used
  )
}
