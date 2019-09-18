.onAttach <- function(...) {
  packageStartupMessage("For execution on a local, multicore CPU with excess RAM we recommend calling")
  packageStartupMessage("options(mc.cores = parallel::detectCores())")
}
