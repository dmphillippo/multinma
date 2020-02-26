
theme_multinma <- function() {
  ggplot2::theme_light() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(colour = "grey70", fill = NA),
      panel.grid.major = ggplot2::element_line(colour = "grey95"),
      panel.grid.minor = ggplot2::element_line(colour = "grey95"),
      strip.background = ggplot2::element_rect(colour = "grey70", fill = "grey90"),
      strip.text = ggplot2::element_text(colour = "black")
    )
}
