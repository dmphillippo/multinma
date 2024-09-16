#' The nma_nodesplit class
#'
#' The `nma_nodesplit` and `nma_nodesplit_df` classes contains the results from
#' running a node-splitting model with the function [nma()].
#'
#' @rdname nma_nodesplit-class
#' @name nma_nodesplit-class
#' @aliases nma_nodesplit nma_nodesplit_df nma_nodesplit_df-class
#'
#' @details Objects of class `nma_nodesplit` inherit from the [stan_nma] class,
#'   and contain the results of fitting a single node-split model. They have one
#'   additional component, `nodesplit`, which gives the comparison that was
#'   node-split as a length 2 vector.
#'
#'   Objects of class `nma_nodesplit_df` are tibble data frames with one row
#'   for each node-split comparison and columns:
#'   \describe{
#'   \item{`trt1`, `trt2`}{Treatments forming the comparison}
#'   \item{`model`}{A list column containing the results of each model as a
#'   `nma_nodesplit` object}
#'   }
#'   Optionally, there will be an additional row for the consistency model if
#'   this was fitted (e.g. by `get_nodesplits(., include_consistency = TRUE)`)
#'   with `trt1` and `trt2` both `NA`.
#'
NULL

#' Print `nma_nodesplit_df` objects
#'
#' @param x A [nma_nodesplit_df] object
#' @param ... Further arguments passed to \code{\link[rstan:print.stanfit]{print.stanfit()}}
#'
#' @seealso The summary method [summary.nma_nodesplit_df()] summarises the
#'   node-splitting results.
#'
#' @return `x` is returned invisibly.
#'
#' @export
print.nma_nodesplit_df <- function(x, ...) {
  n_ns <- nrow(dplyr::filter(x, !is.na(.data$trt1) & !is.na(.data$trt2)))

  cglue("Node-splitting model{if (n_ns > 1) 's' else ''} fitted for {n_ns} comparison{if (n_ns > 1) 's' else ''}.")
  cglue("To summarise these results, use `summary()`.")

  for (i in 1:nrow(x)) {
    cglue("")
    if (n_ns > 1) {
      if (is.na(x$trt1[i]) && is.na(x$trt2[i])) { # Consistency model
        sec_header("Consistency model")
      } else {
        sec_header(glue::glue("Node-split {x$trt2[i]} vs. {x$trt1[i]}"))
      }
      cglue("")
    }
    print(x$model[[i]], ...)
  }

  invisible(x)
}

#' @export
#' @rdname print.nma_nodesplit_df
print.nma_nodesplit <- function(x, ...) {
  NextMethod(...)
}

#' Summarise the results of node-splitting models
#'
#' Posterior summaries of node-splitting models (`nma_nodesplit` and
#' `nma_nodesplit_df` objects) can be produced using the `summary()` method, and
#' plotted using the `plot()` method.
#'
#' @param x,object A `nma_nodesplit` or `nma_nodesplit_df` object
#' @param consistency Optional, a `stan_nma` object for the corresponding fitted
#'   consistency model, to display the network estimates alongside the direct
#'   and indirect estimates. The fitted consistency model present in the
#'   `nma_nodesplit_df` object will be used if this is present (see
#'   [get_nodesplits()]).
#' @param ... Additional arguments passed on to other methods
#' @param probs Numeric vector of specifying quantiles of interest, default
#'   `c(0.025, 0.25, 0.5, 0.75, 0.975)`
#'
#' @details The `plot()` method is a shortcut for `plot(summary(nma_nodesplit))`. For
#'   details of plotting options, see [plot.nodesplit_summary()].
#'
#' @return A [nodesplit_summary] object
#' @export
#'
#' @seealso [plot.nodesplit_summary()]
#'
#' @template ex_smoking_nma_re_nodesplit_example
#' @examples \donttest{
#' # Summarise the node-splitting results
#' summary(smk_fit_RE_nodesplit)
#'
#' # Plot the node-splitting results
#' plot(smk_fit_RE_nodesplit)
#' }
summary.nma_nodesplit_df <- function(object, consistency = NULL, ...,
                                     probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {

  if (any(is.na(object$trt1) & is.na(object$trt2))) { # Consistency model present in object
    consistency <- dplyr::filter(object, is.na(.data$trt1) & is.na(.data$trt2))$model[[1]]

    object <- dplyr::filter(object, !is.na(.data$trt1) | !is.na(.data$trt2))

  } else if (!is.null(consistency)) { # Consistency model provided

    # Check provided consistency object
    if (!inherits(consistency, "stan_nma"))
      abort("`consistency` should be a fitted consistency model (object of class `stan_nma`).")

    if (consistency$consistency != "consistency")
      abort(glue::glue("`consistency` should be a fitted consistency model, not a ",
                       switch(consistency$consistency,
                              nodesplit = "node-splitting",
                              ume = "UME inconsistency"),
                       " model."))

    if (!identical(consistency$stanfit@sim$fnames_oi, setdiff(object$model[[1]]$stanfit@sim$fnames_oi, "omega")))
      abort("The fitted consistency model `consistency` does not match the node-splitting model.")

  }

  out <- object
  out$summary <- purrr::map(out$model, .summarise_ns,
                            cons_mod = consistency, probs = probs)
  out$p_value <- purrr::map_dbl(purrr::map(out$model, as.array, pars = "omega"),
                                ~2 * ifelse(median(.) > 0, mean(. < 0), mean(. > 0)))
  out$dic <- purrr::map(out$model, dic)

  # Drop the model objects from the output
  out$model <- NULL

  class(out) <- c("nodesplit_summary", class(out))

  return(out)
}

#' @export
#' @rdname summary.nma_nodesplit_df
summary.nma_nodesplit <- function(object, consistency = NULL, ...,
                                  probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {

  # Turn nma_nodesplit into nma_nodesplit_df
  df <- tibble::tibble(trt1 = object$nodesplit[1],
                       trt2 = object$nodesplit[2],
                       model = list(object))

  class(df) <- c("nma_nodesplit_df", class(df))

  summary(df, consistency = consistency, ..., probs = probs)
}

#' @export
#' @rdname summary.nma_nodesplit_df
plot.nma_nodesplit <- function(x, consistency = NULL, ...) {

  # Turn nma_nodesplit into nma_nodesplit_df
  df <- tibble::tibble(trt1 = x$nodesplit[1],
                       trt2 = x$nodesplit[2],
                       model = list(x))

  class(df) <- c("nma_nodesplit_df", class(df))

  plot(summary(df, consistency = consistency), ...)
}

#' @export
#' @rdname summary.nma_nodesplit_df
plot.nma_nodesplit_df <- function(x, consistency = NULL, ...) {
  plot(summary(x, consistency = consistency), ...)
}

# Function to summarise direct, indirect, network effects
.summarise_ns <- function(ns_mod, cons_mod = NULL, ..., probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  omega <- as.array(ns_mod, pars = "omega")
  nodesplit <- ns_mod$nodesplit
  ns_text <- paste0(nodesplit[2], " vs. ", nodesplit[1])

  s_dim <- dim(omega)
  s_dim[3] <- 3 # direct, indirect, omega

  s_dimnames <- dimnames(omega)
  s_dimnames[[3]] <- c(paste0("d_dir[", ns_text, "]"),
                       paste0("d_ind[", ns_text, "]"),
                       "omega")

  if (!is.null(cons_mod)) {
    s_dim[3] <- s_dim[3] + 1
    s_dimnames[[3]] <- c(paste0("d_net[", ns_text, "]"),
                         s_dimnames[[3]])
  }
  if (ns_mod$trt_effects == "random") {
    s_dim[3] <- s_dim[3] + 1
    s_dimnames[[3]] <- c(s_dimnames[[3]], "tau")

    if (!is.null(cons_mod)) {
      s_dim[3] <- s_dim[3] + 1
      s_dimnames[[3]] <- c(s_dimnames[[3]], "tau_consistency")
    }
  }

  sims <- array(dim = s_dim, dimnames = s_dimnames)

  if (!is.null(cons_mod)) {
    if (nodesplit[1] == levels(cons_mod$network$treatments)[1]) {
      sims[ , , paste0("d_net[", ns_text, "]")] <-
        as.array(cons_mod, pars = paste0("d[", nodesplit[2], "]"))
    } else {
      sims[ , , paste0("d_net[", ns_text, "]")] <-
        as.array(cons_mod, pars = paste0("d[", nodesplit[2], "]")) -
        as.array(cons_mod, pars = paste0("d[", nodesplit[1], "]"))
    }
  }

  if (nodesplit[1] == levels(ns_mod$network$treatments)[1]) {
    sims[ , , paste0("d_ind[", ns_text, "]")] <-
      as.array(ns_mod, pars = paste0("d[", nodesplit[2], "]"))
  } else {
    sims[ , , paste0("d_ind[", ns_text, "]")] <-
      as.array(ns_mod, pars = paste0("d[", nodesplit[2], "]")) -
      as.array(ns_mod, pars = paste0("d[", nodesplit[1], "]"))
  }

  # Direct estimate = indirect estimate + inconsistency factor omega
  sims[ , , paste0("d_dir[", ns_text, "]")] <- sims[ , , paste0("d_ind[", ns_text, "]"), drop = FALSE] + omega

  sims[ , , "omega"] <- omega

  if(ns_mod$trt_effects == "random") {
    sims[ , , "tau"] <- as.array(ns_mod, pars = "tau")

    if (!is.null(cons_mod)) {
      sims[ , , "tau_consistency"] <- as.array(cons_mod, pars = "tau")
    }
  }

  out <- list(summary = summary_mcmc_array(sims, probs = probs),
              sims = sims)
  class(out) <- "nma_summary"

  return(out)
}
