
nowarn_on_ci <- function(expr, ...) {
  if (isTRUE(as.logical(Sys.getenv("CI")))) {
    suppressWarnings(expr, ...)
  } else {
    expr
  }
}
