.onAttach <- function(...) {
  packageStartupMessage("For execution on a local, multicore CPU with excess RAM we recommend calling")
  packageStartupMessage("options(mc.cores = parallel::detectCores())")

  invisible()
}

.onLoad <- function(...) {
  s3_register("loo::loo", "stan_nma")
  s3_register("loo::waic", "stan_nma")
  s3_register("tidygraph::as_tbl_graph", "nma_data")

  invisible()
}

require_pkg <- function(p, error = TRUE) {
  if (!rlang::is_string(p)) abort("`p` must be a package name as a string")
  if (!rlang::is_bool(error)) abort("`error` must be a logical value TRUE/FALSE")

  if (requireNamespace(p, quietly = TRUE)) {
    invisible(TRUE)
  } else if (error) {
    abort(glue::glue('Install suggested package `{p}` to use this function: install.packages("{p}")'))
  } else {
    invisible(FALSE)
  }
}
