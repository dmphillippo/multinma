.onAttach <- function(...) {
  packageStartupMessage("For execution on a local, multicore CPU with excess RAM we recommend calling")
  packageStartupMessage("options(mc.cores = parallel::detectCores())")

  invisible()
}

.onLoad <- function(...) {
  s3_register("loo::loo", "stan_nma")
  s3_register("loo::waic", "stan_nma")

  invisible()
}
