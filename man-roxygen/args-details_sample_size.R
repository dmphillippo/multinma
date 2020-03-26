#' @details The `sample_size` argument is optional, but when specified:
#'
#'   * Enables centering of predictors (`center = TRUE`) in [nma()] when
#'     a regression model is given
#'   * Enables production of study-specific relative effects, rank probabilities,
#'     etc. for studies in the network when a regression model is given
#'   * Nodes in [plot.nma_data()] may be weighted by sample size
