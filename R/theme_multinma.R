#' Plot theme for multinma plots
#'
#' A simple `ggplot2` theme for plots in the `multinma` package.
#'
#' @param ... Arguments passed to
#'   \code{\link[ggplot2:ggtheme]{ggplot2::theme_light()}}
#'
#' @return A `ggplot2` theme
#' @export
#'
#' @seealso [ggplot2::theme()], \code{\link[ggplot2:theme_get]{ggplot2::theme_set()}}
#'
#' @examples
#' library(ggplot2)
#' theme_set(theme_multinma())
#'
theme_multinma <- function(...) {
  ggplot2::theme_light(...) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "grey70", fill = NA),
      panel.grid.major = ggplot2::element_line(colour = "grey95"),
      panel.grid.minor = ggplot2::element_line(colour = "grey95"),
      strip.background = ggplot2::element_rect(colour = "grey70", fill = "grey90"),
      strip.text = ggplot2::element_text(colour = "black")
    )
}
