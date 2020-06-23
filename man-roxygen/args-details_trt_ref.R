#' @details
#' By default, `trt_ref = NULL` and a network reference treatment will be chosen
#' that attempts to maximise computational efficiency and stability. If an
#' alternative reference treatment is chosen and the model runs slowly or has
#' low effective sample size (ESS) this may be the cause - try letting the
#' default reference treatment be used instead. Regardless of which treatment is
#' used as the network reference at the model fitting stage, results can be
#' transformed afterwards: see the `trt_ref` argument of
#' [relative_effects()] and [predict.stan_nma()].
