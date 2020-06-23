#' Set default values
#'
#' The `.default()` function is used internally to mark certain values as
#' default, so that the user may be notified when default values are being used.
#' For example, choosing a default reference treatment for a network, or using
#' default prior distributions. The function `.is_default()` checks whether an
#' argument/object is set to a default value. Neither of these functions are
#' intended to be called by the user.
#'
#' @param x An object
#'
#' @return For `.default()`, an identical object with additional attribute
#'   `.default`. For `.is_default()`, a logical value (`TRUE` or `FALSE`).
#' @export
#' @rdname default_values
#'
.default <- function(x = list()) {
  attr(x, ".default") <- TRUE
  return(x)
}

#' @export
#' @rdname default_values
.is_default <- function(x) {
  return(isTRUE(attr(x, ".default", exact = TRUE)))
}
