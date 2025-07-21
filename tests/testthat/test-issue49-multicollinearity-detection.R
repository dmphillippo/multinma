# Test file for Issue #49: Enhanced Multicollinearity Detection in add_integration()

test_that("detect_multicollinearity detects high pairwise correlations", {
  # Create a correlation matrix with high correlations
  cor_matrix <- matrix(c(1.0, 0.98, 0.2,
                        0.98, 1.0, 0.1,
                        0.2, 0.1, 1.0), nrow = 3)
  var_names <- c("var1", "var2", "var3")
  
  expect_warning(
    detect_multicollinearity(cor_matrix, var_names),
    "High pairwise correlations detected.*var1 & var2.*r = 0.98"
  )
})

test_that("detect_multicollinearity detects near-singular matrices", {
  # Create a near-singular correlation matrix
  cor_matrix <- matrix(c(1.0, 0.999, 0.999,
                        0.999, 1.0, 0.999,
                        0.999, 0.999, 1.0), nrow = 3)
  var_names <- c("x1", "x2", "x3")
  
  expect_warning(
    detect_multicollinearity(cor_matrix, var_names),
    "Near-singular matrix detected"
  )
})

test_that("detect_multicollinearity detects high condition number", {
  # Create a matrix with high condition number
  cor_matrix <- matrix(c(1.0, 0.95, 0.90,
                        0.95, 1.0, 0.85,
                        0.90, 0.85, 1.0), nrow = 3)
  var_names <- c("cov1", "cov2", "cov3")
  
  expect_warning(
    detect_multicollinearity(cor_matrix, var_names, cond_num_threshold = 10),
    "High condition number detected"
  )
})

test_that("detect_multicollinearity gives positive feedback for well-conditioned matrices", {
  # Create a well-conditioned correlation matrix
  cor_matrix <- matrix(c(1.0, 0.3, 0.2,
                        0.3, 1.0, 0.1,
                        0.2, 0.1, 1.0), nrow = 3)
  var_names <- c("var1", "var2", "var3")
  
  expect_silent(detect_multicollinearity(cor_matrix, var_names))
})

test_that("detect_multicollinearity handles single variable case", {
  cor_matrix <- matrix(1.0)
  var_names <- "single_var"
  
  expect_silent(detect_multicollinearity(cor_matrix, var_names))
})

test_that("detect_multicollinearity validates input parameters", {
  cor_matrix <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  
  # Test invalid cor_matrix
  expect_error(
    detect_multicollinearity("not_a_matrix", c("x", "y")),
    "cor_matrix must be a numeric matrix"
  )
  
  # Test asymmetric matrix
  asym_matrix <- matrix(c(1, 0.5, 0.3, 1), nrow = 2)
  expect_error(
    detect_multicollinearity(asym_matrix, c("x", "y")),
    "cor_matrix must be symmetric"
  )
  
  # Test mismatched var_names length
  expect_error(
    detect_multicollinearity(cor_matrix, "only_one_name"),
    "var_names must be provided and match the number of columns"
  )
})

test_that("multicollinearity detection integrates with add_integration for user-provided correlations", {
  # Create test data
  test_data <- data.frame(
    study = c("Study1", "Study2"),
    y = c(1, 2),
    se = c(0.1, 0.2)
  )
  
  # High correlation matrix
  high_cor <- matrix(c(1.0, 0.98, 0.98, 1.0), nrow = 2)
  
  expect_warning(
    add_integration(test_data, 
                   x1 = distr(qnorm, 0, 1),
                   x2 = distr(qnorm, 0, 1),
                   cor = high_cor),
    "High pairwise correlations detected"
  )
})

test_that("multicollinearity detection integrates with add_integration for IPD-derived correlations", {
  skip_if_not_installed("dplyr")
  
  # Create mock IPD data with high correlations
  ipd_data <- data.frame(
    .study = rep(c("IPD1", "IPD2"), each = 50),
    .trt = sample(c("A", "B"), 100, replace = TRUE),
    x1 = rep(c(1, 2), each = 50) + rnorm(100, 0, 0.01), # Nearly identical
    x2 = rep(c(1, 2), each = 50) + rnorm(100, 0, 0.01), # Nearly identical
    y = rnorm(100)
  )
  
  # Create mock network
  mock_network <- list(
    ipd = ipd_data,
    agd_arm = data.frame(
      .study = "AGD1",
      .trt = "C", 
      y = 1.5,
      se = 0.3
    ),
    agd_contrast = NULL
  )
  class(mock_network) <- "nma_data"
  
  expect_warning(
    add_integration(mock_network,
                   x1 = distr(qnorm, 0, 1),
                   x2 = distr(qnorm, 0, 1)),
    "High pairwise correlations detected"
  )
})

test_that("multicollinearity detection provides helpful remediation advice", {
  cor_matrix <- matrix(c(1.0, 0.98, 0.98, 1.0), nrow = 2)
  var_names <- c("height", "weight")
  
  expect_warning(
    detect_multicollinearity(cor_matrix, var_names),
    "Consider.*removing highly correlated variables.*PCA.*regularization"
  )
})

test_that("multicollinearity detection works with adjusted correlations", {
  # Create test data
  test_data <- data.frame(
    study = c("Study1"),
    y = c(1),
    se = c(0.1)
  )
  
  # Correlation matrix that will cause issues after adjustment
  problematic_cor <- matrix(c(1.0, 0.85, 0.85, 1.0), nrow = 2)
  
  # Should detect issues in adjusted correlations
  expect_warning(
    add_integration(test_data,
                   x1 = distr(qbinom, size = 10, prob = 0.5),
                   x2 = distr(qbinom, size = 10, prob = 0.5),
                   cor = problematic_cor,
                   cor_adjust = "spearman"),
    "multicollinearity|correlation"
  )
}) 