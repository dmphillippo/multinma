#' Bind chains
#'
#' Methods to combine chains from fitted [stan_nma] model objects, or
#' [mcmc_array] 3D MCMC arrays. `cbind.stan_nma()` is an alias for
#' `bind_chains()`.
#'
#' Objects to combine must have the same number of iterations, warmup, and
#' thinning, and all parameter names must match. For `bind_chains()`, this means
#' that both the model specification and input data must be the same.
#'
#' @param ... Multiple fitted [stan_nma] model objects (for `bind_chains()`) or
#' [mcmc_array] arrays (for `cbind.mcmc_array()`). `bind_chains()` also accepts
#' a single list as input, containing the objects to combine.
#'
#' @return A [stan_nma] or [mcmc_array] object.
#' @export
bind_chains <- function(...) {
  m <- list(...)
  if (length(m) == 1 && rlang::is_bare_list(m[[1]])) m <- m[[1]]

  if (!all(purrr::map_lgl(m, inherits, "stan_nma"))) abort("Can only combine stan_nma objects.")

  # check models are compatible
  fnames <- purrr::map(m, ~.$stanfit@sim$fnames_oi)
  if (any(purrr::map_lgl(fnames, ~!identical(., fnames[[1]]))))
    abort(c("Models are not identical.",
            i = "This may be due to different model specifications or different input data."))

  iters <- purrr::map_dfr(m, ~.$stanfit@sim[c("iter", "thin", "warmup")])
  if (any(iters$iter != iters$iter[1])) abort("Number of iterations not equal.")
  if (any(iters$thin != iters$thin[1])) abort("Different thinning values used.")
  if (any(iters$warmup != iters$warmup[1])) abort("Number of warmup iterations not equal.")

  out <- m[[1]]

  out$stanfit <- rstan::sflist2stanfit(purrr::map(m, as.stanfit))

  return(out)
}

#' @rdname bind_chains
#' @export
cbind.stan_nma <- function(...) {
  bind_chains(...)
}

#' @rdname bind_chains
#' @export
cbind.mcmc_array <- function(...) {
  a <- list(...)

  if (!all(purrr::map_lgl(a, inherits, "mcmc_array"))) abort("Can only combine mcmc_array objects.")

  # check dims
  dims <- do.call(rbind, purrr::map(a, dim))
  if (any(dims[, 1] != dims[1, 1])) abort("Number of iterations not equal.")
  if (any(dims[, 3] != dims[1, 3])) abort("Number of parameters not equal.")

  # check names
  nms <- purrr::map(a, ~dimnames(.)[[3]])
  if (any(purrr::map_lgl(nms, ~!identical(., nms[[1]])))) abort("Parameter names do not match.")

  # combine
  odims <- c(dims[1, 1], sum(dims[, 2]), dims[1, 3])
  odimnames <- list(iterations = NULL, chains = NULL, parameters = nms[[1]])

  out <- array(dim = odims, dimnames = odimnames)
  chain <- 0
  for (i in 1:length(a)) {
    out[, chain + (1:dims[i, 2]), ] <- a[[i]]
    chain <- chain + dims[i, 2]
  }

  class(out) <- c("mcmc_array", class(out))
  return(out)
}
