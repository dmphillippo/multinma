#' Target average acceptance probability
#'
#' The Stan control argument `adapt_delta` sets the target average acceptance
#' probability for the No-U-Turn Sampler (NUTS) used by Stan.
#'
#' @details The default value of `adapt_delta` used by [nma()] is 0.8 for fixed
#'   effect models, and 0.95 for random effects models.
#'
#'   You should not need to change `adapt_delta` unless you see a warning
#'   message about divergent transitions. Increasing `adapt_delta` from the
#'   default to a value closer to 1 means that Stan will use a smaller step
#'   size, making sampling slower but more robust, and resulting in fewer
#'   divergent transitions.
#'
#'   For more details see the Stan documentation available from
#'   \url{https://mc-stan.org/users/documentation/}.
#'
#' @name adapt_delta
#'
NULL
