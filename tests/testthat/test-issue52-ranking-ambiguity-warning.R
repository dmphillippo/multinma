# Test file for Issue #52: TTE Parametric Distribution Ranking Warnings

# Mock stan_nma object helper function
create_mock_stan_nma <- function(likelihood = "normal", aux_by = NULL, aux_regression = NULL) {
  # Mock network structure
  mock_network <- list(
    treatments = factor(c("A", "B", "C"))
  )
  
  # Mock stan_nma object
  stan_nma_mock <- list(
    likelihood = likelihood,
    consistency = "consistency",
    network = mock_network,
    aux_by = aux_by,
    aux_regression = aux_regression
  )
  class(stan_nma_mock) <- "stan_nma"
  
  return(stan_nma_mock)
}

# Mock relative_effects function
with_mock_relative_effects <- function(test_expr) {
  mockery::stub(posterior_ranks, "relative_effects", function(x, newdata, study, all_contrasts, summary) {
    # Return a mock relative effects object
    list(
      sim = array(rnorm(3 * 4 * 2), dim = c(3, 4, 2), 
                  dimnames = list(iterations = NULL, chains = NULL, 
                                 parameters = c("d[B]", "d[C]"))),
      studies = NULL
    )
  })
  
  test_expr
}

test_that("posterior_ranks warns for TTE parametric models with aux_by treatment effects", {
  # Test various TTE likelihoods with aux_by including .trt
  tte_likelihoods <- c("weibull", "gompertz", "weibull-aft", "lognormal", 
                      "loglogistic", "gamma", "gengamma")
  
  for (likelihood in tte_likelihoods) {
    mock_obj <- create_mock_stan_nma(likelihood = likelihood, aux_by = c("age", ".trt"))
    
    expect_warning(
      with_mock_relative_effects({
        posterior_ranks(mock_obj)
      }),
      "Treatment rankings may be ambiguous.*time-to-event parametric models.*auxiliary treatment effects"
    )
  }
})

test_that("posterior_ranks warns for TTE parametric models with aux_regression treatment effects", {
  # Create mock aux_regression formula with .trt terms
  mock_formula <- y ~ age + .trt
  attr(mock_formula, "factors") <- matrix(c(0, 1, 0, 0, 0, 1), nrow = 3, 
                                         dimnames = list(c("y", "age", ".trt"), c("age", ".trt")))
  
  mock_obj <- create_mock_stan_nma(likelihood = "weibull", aux_regression = mock_formula)
  
  expect_warning(
    with_mock_relative_effects({
      posterior_ranks(mock_obj)
    }),
    "Treatment rankings may be ambiguous.*time-to-event parametric models.*auxiliary treatment effects"
  )
})

test_that("posterior_ranks does not warn for non-TTE likelihoods with treatment effects", {
  non_tte_likelihoods <- c("normal", "binomial", "poisson", "ordered_multinomial")
  
  for (likelihood in non_tte_likelihoods) {
    mock_obj <- create_mock_stan_nma(likelihood = likelihood, aux_by = c("age", ".trt"))
    
    expect_silent(
      with_mock_relative_effects({
        suppressMessages(posterior_ranks(mock_obj))
      })
    )
  }
})

test_that("posterior_ranks does not warn for TTE likelihoods without treatment effects", {
  tte_likelihoods <- c("weibull", "gompertz", "weibull-aft", "lognormal", 
                      "loglogistic", "gamma", "gengamma")
  
  for (likelihood in tte_likelihoods) {
    # Test without aux_by
    mock_obj1 <- create_mock_stan_nma(likelihood = likelihood, aux_by = c("age", "sex"))
    expect_silent(
      with_mock_relative_effects({
        suppressMessages(posterior_ranks(mock_obj1))
      })
    )
    
    # Test without aux_regression
    mock_obj2 <- create_mock_stan_nma(likelihood = likelihood, aux_regression = NULL)
    expect_silent(
      with_mock_relative_effects({
        suppressMessages(posterior_ranks(mock_obj2))
      })
    )
    
    # Test with aux_regression but no .trt terms
    mock_formula <- y ~ age + sex
    attr(mock_formula, "factors") <- matrix(c(0, 1, 0, 0, 0, 1), nrow = 3, 
                                           dimnames = list(c("y", "age", "sex"), c("age", "sex")))
    mock_obj3 <- create_mock_stan_nma(likelihood = likelihood, aux_regression = mock_formula)
    expect_silent(
      with_mock_relative_effects({
        suppressMessages(posterior_ranks(mock_obj3))
      })
    )
  }
})

test_that("posterior_ranks warning includes helpful guidance", {
  mock_obj <- create_mock_stan_nma(likelihood = "weibull", aux_by = c("age", ".trt"))
  
  expect_warning(
    with_mock_relative_effects({
      posterior_ranks(mock_obj)
    }),
    "Rankings are based on main treatment effects.*but auxiliary treatment effects.*also influence.*Consider examining both sets of effects"
  )
})

test_that("posterior_ranks handles edge cases correctly", {
  # Test with NULL aux_by and aux_regression
  mock_obj1 <- create_mock_stan_nma(likelihood = "weibull", aux_by = NULL, aux_regression = NULL)
  expect_silent(
    with_mock_relative_effects({
      suppressMessages(posterior_ranks(mock_obj1))
    })
  )
  
  # Test with empty aux_by
  mock_obj2 <- create_mock_stan_nma(likelihood = "weibull", aux_by = character(0))
  expect_silent(
    with_mock_relative_effects({
      suppressMessages(posterior_ranks(mock_obj2))
    })
  )
  
  # Test with .trt not at the beginning of aux_regression terms
  mock_formula <- y ~ age + trt_class + .trt_interaction
  attr(mock_formula, "factors") <- matrix(c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), nrow = 4, 
                                         dimnames = list(c("y", "age", "trt_class", ".trt_interaction"), 
                                                        c("age", "trt_class", ".trt_interaction")))
  mock_obj3 <- create_mock_stan_nma(likelihood = "weibull", aux_regression = mock_formula)
  expect_silent(
    with_mock_relative_effects({
      suppressMessages(posterior_ranks(mock_obj3))
    })
  )
})

test_that("posterior_ranks warning detection logic is robust", {
  # Test case where aux_by contains .trt among other variables
  mock_obj1 <- create_mock_stan_nma(likelihood = "gamma", aux_by = c("baseline_risk", ".trt", "study_year"))
  expect_warning(
    with_mock_relative_effects({
      posterior_ranks(mock_obj1)
    }),
    "Treatment rankings may be ambiguous"
  )
  
  # Test case with complex aux_regression including .trt terms
  mock_formula <- y ~ baseline + .trt + baseline:.trt
  attr(mock_formula, "factors") <- matrix(c(0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1), nrow = 4, 
                                         dimnames = list(c("y", "baseline", ".trt", "baseline:.trt"), 
                                                        c("baseline", ".trt", "baseline:.trt", ".trt")))
  mock_obj2 <- create_mock_stan_nma(likelihood = "gengamma", aux_regression = mock_formula)
  expect_warning(
    with_mock_relative_effects({
      posterior_ranks(mock_obj2)
    }),
    "Treatment rankings may be ambiguous"
  )
}) 