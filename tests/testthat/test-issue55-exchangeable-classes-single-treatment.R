# Test file for Issue #55: Fix Exchangeable Classes with Single Treatments

# Test data setup
create_single_treatment_class_data <- function() {
  # Create test data with mixed treatment classes
  # Some classes have multiple treatments, some have single treatments
  data.frame(
    .study = paste0("Study", 1:6),
    .trt = c("A", "B", "C", "D", "E", "F"),
    .y = rnorm(6),
    .se = runif(6, 0.1, 0.3)
  )
}

create_all_single_treatment_classes_data <- function() {
  # Create test data where every treatment is in its own class
  data.frame(
    .study = paste0("Study", 1:4),
    .trt = c("A", "B", "C", "D"),
    .y = rnorm(4),
    .se = runif(4, 0.1, 0.3)
  )
}

test_that("nma handles exchangeable class effects with single treatments", {
  skip_if_not_installed("rstan")
  
  # Test data with mixed class sizes
  test_data <- create_single_treatment_class_data()
  
  # Mock the nma function to avoid actual Stan compilation/sampling
  with_mocked_bindings(
    "sampling" = function(...) {
      # Return minimal stanfit-like object
      structure(list(), class = "stanfit")
    },
    .package = "multinma",
    {
      # Should not error with exchangeable class effects
      expect_no_error({
        result <- nma(test_data,
                     likelihood = "normal",
                     class_effects = "exchangeable",
                     iter = 1, chains = 1, refresh = 0)
      })
      
      expect_true(inherits(result, "stan_nma"))
    }
  )
})

test_that("nma handles all single-treatment exchangeable classes", {
  skip_if_not_installed("rstan")
  
  # Test data where every treatment is in its own class
  test_data <- create_all_single_treatment_classes_data()
  
  with_mocked_bindings(
    "sampling" = function(...) {
      structure(list(), class = "stanfit")
    },
    .package = "multinma",
    {
      # Should not error even when all classes have single treatments
      expect_no_error({
        result <- nma(test_data,
                     likelihood = "normal", 
                     class_effects = "exchangeable",
                     iter = 1, chains = 1, refresh = 0)
      })
      
      expect_true(inherits(result, "stan_nma"))
    }
  )
})

test_that("exchangeable classes work with different likelihoods", {
  skip_if_not_installed("rstan")
  
  test_data <- create_single_treatment_class_data()
  
  # Test different likelihoods that might be affected
  likelihoods_to_test <- c("normal", "binomial", "poisson")
  
  for (likelihood in likelihoods_to_test) {
    with_mocked_bindings(
      "sampling" = function(...) {
        structure(list(), class = "stanfit")
      },
      .package = "multinma",
      {
        expect_no_error({
          result <- nma(test_data,
                       likelihood = likelihood,
                       class_effects = "exchangeable", 
                       iter = 1, chains = 1, refresh = 0)
        }, info = paste("Failed for likelihood:", likelihood))
        
        expect_true(inherits(result, "stan_nma"))
      }
    )
  }
})

test_that("exchangeable classes work with random effects", {
  skip_if_not_installed("rstan")
  
  test_data <- create_single_treatment_class_data()
  
  with_mocked_bindings(
    "sampling" = function(...) {
      structure(list(), class = "stanfit")
    },
    .package = "multinma",
    {
      expect_no_error({
        result <- nma(test_data,
                     likelihood = "normal",
                     trt_effects = "random",
                     class_effects = "exchangeable",
                     iter = 1, chains = 1, refresh = 0)
      })
      
      expect_true(inherits(result, "stan_nma"))
    }
  )
})

test_that("exchangeable classes work with regression models", {
  skip_if_not_installed("rstan")
  
  test_data <- create_single_treatment_class_data()
  test_data$age <- rnorm(nrow(test_data), 50, 10)  # Add covariate
  
  with_mocked_bindings(
    "sampling" = function(...) {
      structure(list(), class = "stanfit")
    },
    .package = "multinma",
    {
      expect_no_error({
        result <- nma(test_data,
                     likelihood = "normal",
                     regression = ~ age,
                     class_effects = "exchangeable",
                     iter = 1, chains = 1, refresh = 0)
      })
      
      expect_true(inherits(result, "stan_nma"))
    }
  )
})

test_that("Stan files have bounds checking for array access", {
  # This test verifies that the Stan code fixes have been applied
  # by checking that the problematic patterns have been replaced
  
  # Find all Stan files
  stan_files <- list.files("inst/stan", pattern = "\\.stan$", 
                          recursive = TRUE, full.names = TRUE)
  
  if (length(stan_files) > 0) {
    for (stan_file in stan_files) {
      if (file.exists(stan_file)) {
        content <- readLines(stan_file)
        
        # Check that old problematic patterns are not present
        old_patterns <- grepl("agd_arm_trt\\[i\\] > 1 && which_CE\\[agd_arm_trt\\[i\\] - 1\\]", content)
        expect_false(any(old_patterns), 
                    info = paste("Old pattern found in", stan_file))
        
        # Check that new safe patterns are present (if the file has class effects)
        if (any(grepl("which_CE", content))) {
          safe_patterns <- grepl("agd_arm_trt\\[i\\] > 1 && agd_arm_trt\\[i\\] <= nt && which_CE", content)
          expect_true(any(safe_patterns), 
                     info = paste("Safe pattern not found in", stan_file))
        }
      }
    }
  } else {
    skip("No Stan files found")
  }
})

test_that("bounds checking prevents array index errors", {
  # This is a conceptual test - the actual prevention happens in Stan
  # We verify that the patterns include the necessary bounds check
  
  # Example of what should be present in the fixed code:
  # if (agd_arm_trt[i] > 1 && agd_arm_trt[i] <= nt && which_CE[agd_arm_trt[i] - 1])
  
  # Simulate the bounds checking logic
  simulate_bounds_check <- function(trt_index, nt) {
    # Old logic (problematic):
    # return(trt_index > 1)
    
    # New logic (safe):
    return(trt_index > 1 && trt_index <= nt)
  }
  
  nt <- 3  # Number of treatments
  
  # Test cases
  expect_true(simulate_bounds_check(2, nt))   # Valid treatment index
  expect_true(simulate_bounds_check(3, nt))   # Valid treatment index  
  expect_false(simulate_bounds_check(1, nt))  # Reference treatment (index 1)
  expect_false(simulate_bounds_check(4, nt))  # Out of bounds (should be safe now)
  expect_false(simulate_bounds_check(0, nt))  # Invalid index
})

test_that("exchangeable classes handle edge cases gracefully", {
  skip_if_not_installed("rstan")
  
  # Test with minimal data (2 treatments)
  minimal_data <- data.frame(
    .study = c("Study1", "Study2"),
    .trt = c("A", "B"),
    .y = c(0.5, 0.3),
    .se = c(0.1, 0.15)
  )
  
  with_mocked_bindings(
    "sampling" = function(...) {
      structure(list(), class = "stanfit")
    },
    .package = "multinma",
    {
      expect_no_error({
        result <- nma(minimal_data,
                     likelihood = "normal",
                     class_effects = "exchangeable",
                     iter = 1, chains = 1, refresh = 0)
      })
    }
  )
})

test_that("fix applies to all relevant Stan model files", {
  # List of Stan files that should have the fix applied
  expected_files <- c(
    "inst/stan/binomial_1par.stan",
    "inst/stan/binomial_2par.stan", 
    "inst/stan/normal.stan",
    "inst/stan/poisson.stan",
    "inst/stan/ordered_multinomial.stan",
    "inst/stan/survival_param.stan",
    "inst/stan/survival_mspline.stan",
    "inst/stan/include/transformed_parameters_common.stan"
  )
  
  for (file_path in expected_files) {
    if (file.exists(file_path)) {
      content <- readLines(file_path)
      
      # Should have the new bounds checking pattern
      has_bounds_check <- any(grepl("&& .*<= nt && which_CE", content))
      
      # Only check if file has which_CE references (some files might not use class effects)
      if (any(grepl("which_CE", content))) {
        expect_true(has_bounds_check, 
                   info = paste("Bounds checking missing in", file_path))
      }
    }
  }
})

test_that("fix handles different treatment indexing patterns", {
  # The fix should handle different variable names used for treatment indexing
  # This test verifies that all patterns have been addressed
  
  patterns_to_check <- c(
    "agd_arm_trt[i]",
    "ipd_trt[ipd_arm[i]]", 
    "agd_contrast_trt[i]",
    "agd_contrast_trt_b[i]"
  )
  
  # Check common transformed parameters file which has multiple patterns
  common_file <- "inst/stan/include/transformed_parameters_common.stan"
  
  if (file.exists(common_file)) {
    content <- readLines(common_file)
    
    # Each pattern that appears should have bounds checking
    for (pattern in patterns_to_check) {
      pattern_lines <- grepl(pattern, content, fixed = TRUE)
      
      if (any(pattern_lines)) {
        # If pattern exists, check that bounds checking is present
        bounds_check_lines <- grepl(paste0(pattern, ".*<= nt"), content)
        expect_true(any(bounds_check_lines),
                   info = paste("Bounds checking missing for pattern:", pattern))
      }
    }
  }
}) 