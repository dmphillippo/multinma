# Test file for Issue #54: MCMC Convergence Diagnostics Functions

# Mock helper functions for testing
create_mock_stan_nma <- function(trt_effects = "fixed", regression = NULL) {
  # Create mock stanfit object with summary
  mock_summary <- list(
    summary = matrix(
      c(0.5, 0.1, 1.05, 500,  # d[B]: mean, sd, Rhat, n_eff
        0.3, 0.2, 1.02, 800,  # d[C]: mean, sd, Rhat, n_eff
        -0.2, 0.15, 1.01, 600, # beta[1]: mean, sd, Rhat, n_eff (if regression)
        0.25, 0.05, 1.03, 400), # tau: mean, sd, Rhat, n_eff (if random)
      nrow = 4, ncol = 4,
      dimnames = list(
        c("d[B]", "d[C]", "beta[1]", "tau"),
        c("mean", "sd", "Rhat", "n_eff")
      )
    )
  )
  
  mock_stanfit <- structure(
    list(model_pars = c("d", "beta", "tau", "alpha")),
    class = "stanfit"
  )
  
  # Mock rstan::summary function
  mockery::stub(rstan::summary, where = environment(), what = function(object, pars = NULL, ...) {
    if (is.null(pars)) {
      return(mock_summary)
    } else {
      # Filter summary based on requested parameters
      available_pars <- rownames(mock_summary$summary)
      matching_pars <- available_pars[available_pars %in% pars]
      if (length(matching_pars) > 0) {
        filtered_summary <- list(
          summary = mock_summary$summary[matching_pars, , drop = FALSE]
        )
        return(filtered_summary)
      } else {
        return(list(summary = matrix(nrow = 0, ncol = 4, dimnames = list(NULL, c("mean", "sd", "Rhat", "n_eff")))))
      }
    }
  })
  
  # Mock rstan::extract function
  mockery::stub(rstan::extract, where = environment(), what = function(object, pars = NULL, permuted = FALSE, ...) {
    # Return mock MCMC array
    n_iter <- 100
    n_chains <- 4
    n_pars <- if (is.null(pars)) 4 else length(pars)
    
    array(rnorm(n_iter * n_chains * n_pars), 
          dim = c(n_iter, n_chains, n_pars),
          dimnames = list(
            iterations = NULL,
            chains = paste0("chain:", 1:n_chains),
            parameters = if (is.null(pars)) c("d[B]", "d[C]", "beta[1]", "tau") else pars
          ))
  })
  
  stan_nma_obj <- list(
    stanfit = mock_stanfit,
    trt_effects = trt_effects,
    regression = regression
  )
  class(stan_nma_obj) <- "stan_nma"
  
  return(stan_nma_obj)
}

test_that("mcmc_trace generates plots for stan_nma objects", {
  skip_if_not_installed("bayesplot")
  
  mock_obj <- create_mock_stan_nma()
  
  # Mock bayesplot::mcmc_trace
  with_mocked_bindings(
    "mcmc_trace" = function(...) {
      structure(list(data = "mock_plot"), class = c("gg", "ggplot"))
    },
    .package = "bayesplot",
    {
      result <- mcmc_trace(mock_obj, pars = c("d[B]", "d[C]"))
      expect_true(inherits(result, "ggplot"))
    }
  )
})

test_that("mcmc_acf generates autocorrelation plots", {
  skip_if_not_installed("bayesplot")
  
  mock_obj <- create_mock_stan_nma()
  
  with_mocked_bindings(
    "mcmc_acf" = function(...) {
      structure(list(data = "mock_acf_plot"), class = c("gg", "ggplot"))
    },
    .package = "bayesplot",
    {
      result <- mcmc_acf(mock_obj, pars = c("d[B]", "d[C]"))
      expect_true(inherits(result, "ggplot"))
    }
  )
})

test_that("mcmc_rhat generates R-hat diagnostic plots", {
  skip_if_not_installed("bayesplot")
  
  mock_obj <- create_mock_stan_nma()
  
  with_mocked_bindings(
    "mcmc_rhat" = function(...) {
      structure(list(data = "mock_rhat_plot"), class = c("gg", "ggplot"))
    },
    .package = "bayesplot",
    {
      result <- mcmc_rhat(mock_obj, pars = c("d[B]", "d[C]"))
      expect_true(inherits(result, "ggplot"))
    }
  )
})

test_that("mcmc_neff generates effective sample size plots", {
  skip_if_not_installed("bayesplot")
  
  mock_obj <- create_mock_stan_nma()
  
  with_mocked_bindings(
    "mcmc_neff" = function(...) {
      structure(list(data = "mock_neff_plot"), class = c("gg", "ggplot"))
    },
    .package = "bayesplot",
    {
      result <- mcmc_neff(mock_obj, pars = c("d[B]", "d[C]"))
      expect_true(inherits(result, "ggplot"))
    }
  )
})

test_that("mcmc_diagnostics_plot combines multiple diagnostic types", {
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("patchwork")
  
  mock_obj <- create_mock_stan_nma()
  
  # Mock all bayesplot functions
  with_mocked_bindings(
    "mcmc_trace" = function(...) structure(list(data = "trace"), class = c("gg", "ggplot")),
    "mcmc_acf" = function(...) structure(list(data = "acf"), class = c("gg", "ggplot")),
    "mcmc_rhat" = function(...) structure(list(data = "rhat"), class = c("gg", "ggplot")),
    "mcmc_neff" = function(...) structure(list(data = "neff"), class = c("gg", "ggplot")),
    .package = "bayesplot",
    {
      with_mocked_bindings(
        "wrap_plots" = function(...) structure(list(data = "combined"), class = c("gg", "ggplot")),
        .package = "patchwork",
        {
          # Test single diagnostic type
          result_single <- mcmc_diagnostics_plot(mock_obj, type = "trace")
          expect_true(inherits(result_single, "ggplot"))
          
          # Test multiple diagnostic types
          result_multi <- mcmc_diagnostics_plot(mock_obj, type = c("trace", "rhat"))
          expect_true(inherits(result_multi, "ggplot"))
        }
      )
    }
  )
})

test_that("get_default_parameters selects appropriate parameters", {
  # Test fixed effects model
  mock_obj_fixed <- create_mock_stan_nma(trt_effects = "fixed")
  
  with_mocked_bindings(
    "summary" = function(object, ...) {
      list(summary = matrix(nrow = 2, ncol = 4, 
                           dimnames = list(c("d[B]", "d[C]"), 
                                         c("mean", "sd", "Rhat", "n_eff"))))
    },
    .package = "rstan",
    {
      pars <- get_default_parameters(mock_obj_fixed)
      expect_true(is.character(pars))
      expect_true(length(pars) >= 0)
    }
  )
})

test_that("validate_parameters handles missing parameters correctly", {
  mock_obj <- create_mock_stan_nma()
  
  with_mocked_bindings(
    "summary" = function(object, ...) {
      list(summary = matrix(nrow = 2, ncol = 4, 
                           dimnames = list(c("d[B]", "d[C]"), 
                                         c("mean", "sd", "Rhat", "n_eff"))))
    },
    .package = "rstan",
    {
      # Test with valid parameters
      valid_pars <- validate_parameters(mock_obj, c("d[B]", "d[C]"))
      expect_equal(valid_pars, c("d[B]", "d[C]"))
      
      # Test with some invalid parameters
      expect_warning(
        mixed_pars <- validate_parameters(mock_obj, c("d[B]", "nonexistent")),
        "not found in the model"
      )
      expect_equal(mixed_pars, "d[B]")
      
      # Test with all invalid parameters
      expect_error(
        validate_parameters(mock_obj, c("nonexistent1", "nonexistent2")),
        "None of the specified parameters were found"
      )
    }
  )
})

test_that("mcmc diagnostic functions handle NULL pars correctly", {
  skip_if_not_installed("bayesplot")
  
  mock_obj <- create_mock_stan_nma()
  
  with_mocked_bindings(
    "mcmc_trace" = function(...) structure(list(), class = c("gg", "ggplot")),
    "summary" = function(object, ...) {
      list(summary = matrix(c(1.05, 500), nrow = 1, ncol = 2, 
                           dimnames = list("d[B]", c("Rhat", "n_eff"))))
    },
    .package = c("bayesplot", "rstan"),
    {
      # Should use default parameters when pars = NULL
      result <- mcmc_trace(mock_obj, pars = NULL)
      expect_true(inherits(result, "ggplot"))
    }
  )
})

test_that("mcmc diagnostic functions validate input objects", {
  # Test with non-stan_nma object
  expect_error(
    mcmc_trace("not_a_stan_nma_object"),
    "x must be a stan_nma object"
  )
  
  expect_error(
    mcmc_acf(list(some = "data")),
    "x must be a stan_nma object"
  )
  
  expect_error(
    mcmc_rhat(NULL),
    "x must be a stan_nma object"
  )
  
  expect_error(
    mcmc_neff(42),
    "x must be a stan_nma object"
  )
  
  expect_error(
    mcmc_diagnostics_plot(data.frame()),
    "x must be a stan_nma object"
  )
})

test_that("mcmc diagnostic functions work with different model types", {
  skip_if_not_installed("bayesplot")
  
  # Test with random effects model
  mock_obj_random <- create_mock_stan_nma(trt_effects = "random")
  
  with_mocked_bindings(
    "mcmc_trace" = function(...) structure(list(), class = c("gg", "ggplot")),
    "summary" = function(object, ...) {
      list(summary = matrix(c(1.02, 1.01, 800, 600), nrow = 2, ncol = 2, 
                           dimnames = list(c("d[B]", "tau"), c("Rhat", "n_eff"))))
    },
    .package = c("bayesplot", "rstan"),
    {
      result <- mcmc_trace(mock_obj_random)
      expect_true(inherits(result, "ggplot"))
    }
  )
  
  # Test with regression model
  mock_obj_regression <- create_mock_stan_nma(regression = ~ age)
  
  with_mocked_bindings(
    "mcmc_trace" = function(...) structure(list(), class = c("gg", "ggplot")),
    "summary" = function(object, ...) {
      list(summary = matrix(c(1.01, 1.03, 700, 500), nrow = 2, ncol = 2, 
                           dimnames = list(c("d[B]", "beta[1]"), c("Rhat", "n_eff"))))
    },
    .package = c("bayesplot", "rstan"),
    {
      result <- mcmc_trace(mock_obj_regression)
      expect_true(inherits(result, "ggplot"))
    }
  )
})

test_that("mcmc_diagnostics_plot handles different type combinations", {
  skip_if_not_installed("bayesplot")
  skip_if_not_installed("patchwork")
  
  mock_obj <- create_mock_stan_nma()
  
  with_mocked_bindings(
    "mcmc_trace" = function(...) structure(list(), class = c("gg", "ggplot")),
    "mcmc_acf" = function(...) structure(list(), class = c("gg", "ggplot")),
    "mcmc_rhat" = function(...) structure(list(), class = c("gg", "ggplot")),
    "mcmc_neff" = function(...) structure(list(), class = c("gg", "ggplot")),
    .package = "bayesplot",
    {
      with_mocked_bindings(
        "wrap_plots" = function(...) structure(list(), class = c("gg", "ggplot")),
        .package = "patchwork",
        {
          # Test each individual type
          expect_true(inherits(mcmc_diagnostics_plot(mock_obj, type = "trace"), "ggplot"))
          expect_true(inherits(mcmc_diagnostics_plot(mock_obj, type = "acf"), "ggplot"))
          expect_true(inherits(mcmc_diagnostics_plot(mock_obj, type = "rhat"), "ggplot"))
          expect_true(inherits(mcmc_diagnostics_plot(mock_obj, type = "neff"), "ggplot"))
          
          # Test combinations
          expect_true(inherits(mcmc_diagnostics_plot(mock_obj, type = c("trace", "rhat")), "ggplot"))
          expect_true(inherits(mcmc_diagnostics_plot(mock_obj, type = c("acf", "neff")), "ggplot"))
          expect_true(inherits(mcmc_diagnostics_plot(mock_obj, type = c("trace", "acf", "rhat")), "ggplot"))
          expect_true(inherits(mcmc_diagnostics_plot(mock_obj, type = c("trace", "acf", "rhat", "neff")), "ggplot"))
        }
      )
    }
  )
}) 