
if (requireNamespace("pkgdown", quietly = TRUE) && pkgdown::in_pkgdown()) { # Setup for pkgdown articles
  knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.align = "center"
  )
  options(width = 95, rstan_refresh = 0)

} else { # Setup for vignettes
  knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    dev = "ragg_png",
    dpi = 96,
    fig.width = 6,
    fig.height = 4.5,
    fig.align = "center",
    fig.cap = "",
    pngquant = '--speed 4 --skip-if-larger'
  )
  options(width = 100, rstan_refresh = 0)

  # Compress PNG images - needs pngquant on path
  knitr::knit_hooks$set(pngquant = knitr::hook_pngquant)

  # Write vignette details to asis file for precompiled vignettes
  if (!rlang::has_name(rmarkdown::metadata, "vignette")) {
    .asis <- paste0(gsub("\\.Rmd(\\.orig)?$", "", knitr::current_input()), ".html.asis")
    file.create(.asis)
    vig_title <- rmarkdown::metadata$title
    if (is.null(vig_title)) vig_title <- rmarkdown::yaml_front_matter(knitr::current_input())$title
    cat(
      paste0("%\\VignetteIndexEntry{", vig_title, "}"),
      "%\\VignetteEngine{R.rsp::asis}",
      "%\\VignetteEncoding{UTF-8}",
      sep = "\n",
      file = .asis)
  }
}

nowarn_on_ci <- function(expr, ...) {
  if (isTRUE(as.logical(Sys.getenv("CI")))) {
    suppressWarnings(expr, ...)
  } else {
    expr
  }
}
