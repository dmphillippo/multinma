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
#' @return For `.default()`, an identical object with additional class
#'   `.default`. For `.is_default()`, a logical value (`TRUE` or `FALSE`).
#' @export
#' @rdname default_values
#'
.default <- function(x = list()) {
  x_call <- rlang::enquo(x)
  return(structure(x, class = c(".default", class(x)), call = x_call))
}

#' @export
#' @rdname default_values
.is_default <- function(x) {
  return(inherits(x, ".default"))
}

get_default_call <- function(x) {
  return(rlang::as_label(attr(x, "call")))
}
