# print.R
#' Print a synpdx_fit object
#'
#' @param x A synpdx_fit object.
#' @param ... Unused.
#' @export
print.synpdx_fit <- function(x, ...) {
  cat("synpdx_fit object\n")
  mdl <- x$selected_model %||% x$best_model_name %||% NA_character_
  if (!is.na(mdl)) cat(" model:", mdl, "\n")
  if (!is.null(x$saemix)) {
    cat(" saemix summary:\n")
    out <- try(utils::capture.output(show(x$saemix)), silent = TRUE)
    if (!inherits(out, "try-error")) cat(paste(out, collapse = "\n"), "\n")
  }
  invisible(x)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
