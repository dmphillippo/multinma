
test_that("`inclusive` argument is validated", {
  m <- "`inclusive` must be a logical value"

  expect_error(multi(a = 1:4, b = 1:4, inclusive = "a"), m)
  expect_error(multi(a = 1:4, b = 1:4, inclusive = 1), m)
  expect_error(multi(a = 1:4, b = 1:4, inclusive = list()), m)
  expect_error(multi(a = 1:4, b = 1:4, inclusive = NA), m)
  expect_error(multi(a = 1:4, b = 1:4, inclusive = NULL), m)

  expect_warning(multi(1, 2 , inclusive = TRUE, type = "competing"),
                 "Ignoring inclusive = TRUE")
})

test_that("`type` argument is validated", {
  m <- "`type` must be"

  expect_error(multi(a = 1:4, b = 1:4, type = "a"), m)
  expect_error(multi(a = 1:4, b = 1:4, type = 1), m)
  expect_error(multi(a = 1:4, b = 1:4, type = list()), m)
  expect_error(multi(a = 1:4, b = 1:4, type = NA), m)
  expect_error(multi(a = 1:4, b = 1:4, type = NULL), m)
})

test_that("basic validation of `...`", {
  expect_error(multi(), "At least 2 outcomes")
  expect_error(multi(1:3), "At least 2 outcomes")
  expect_error(multi(a = 1:4, a = 1:4), "Duplicate outcome category labels")
  a <- 1:4
  expect_error(multi(a, a), "Duplicate outcome category labels")
  expect_error(multi(1:2, 1:4), "must be the same length")
  expect_error(multi(a = 1, b = Inf), "cannot be Inf")
  expect_error(multi(a = 1, b = NaN), "cannot be NaN")
  expect_error(multi(a = 1, b = "a"), "must be numeric")
  expect_error(multi(a = 1, b = 1.5), "must be integer")
  expect_error(multi(a = 1, b = -1), "must be non-negative")
})

df_inclusive <- tibble::tribble(~a, ~b, ~c,
                                1, 1, 1,
                                5, 4, 1,
                                5, 2, 2,
                                10, 5, 0,
                                5, 5, 0,    # y > -Inf, y > b, y > c
                                7, NA, 6,   # y > -Inf, y > c
                                10, 4, NA)  # y > -Inf, y > b

# internally stored as exclusive counts
df_exclusive <- tibble::tribble(~a, ~b, ~c,
                                0, 0, 1,
                                1, 3, 1,
                                3, 0, 2,
                                5, 5, 0,
                                0, 5, 0,    # -Inf < y < b, b < y < c, y > c
                                1, NA, 6,   # -Inf < y < c, y > c
                                6, 4, NA)   # -Inf < y < b, y > b

test_that("ordered inclusive outcome validation", {
  m <- "must be decreasing or constant across increasing categories"

  expect_error(multi(1, 2, 3, inclusive = TRUE, type = "ordered"), paste0(m, ".+row 1"))
  expect_error(multi(5, 4, 6, inclusive = TRUE, type = "ordered"), paste0(m, ".+row 1"))
  expect_error(multi(5, NA, 6, inclusive = TRUE, type = "ordered"), paste0(m, ".+row 1"))
  expect_error(with(df_inclusive[-7, ], multi(c, b, a, inclusive = TRUE, type = "ordered")), paste0(m, ".+rows 2, 3, 4, 5 and 6"))

  expect_error(multi(a = NA, b = 10, c = 5, inclusive = TRUE, type = "ordered"),
               "cannot be missing in the lowest category.+row 1")

  expect_error(multi(a = 10, b = NA, inclusive = TRUE, type = "ordered"),
               "must be present for at least 2 categories.+row 1")

  o <- with(df_inclusive, multi(a, b, c, inclusive = TRUE, type = "ordered"))
  expect_s3_class(o, "multi_ordered")
  expect_equivalent(unclass(o), as.matrix(df_exclusive))
  expect_equal(colnames(o), colnames(df_inclusive))

  o2 <- with(df_inclusive, multi(A = a, B = b, C = c, inclusive = TRUE, type = "ordered"))
  expect_equal(colnames(o2), c("A", "B", "C"))
})

test_that("exclusive outcome validation", {
  expect_error(multi(a = NA, b = 10, c = 5, inclusive = FALSE, type = "ordered"),
               "cannot be missing in the lowest category.+row 1")

  expect_error(multi(a = 10, b = NA, inclusive = FALSE, type = "ordered"),
               "must be present for at least 2 categories.+row 1")
  expect_error(multi(a = 10, b = NA, inclusive = FALSE, type = "competing"),
               "must be present for at least 2 categories.+row 1")

  o <- with(df_exclusive, multi(a, b, c, inclusive = FALSE, type = "ordered"))
  expect_s3_class(o, "multi_ordered")
  expect_equivalent(unclass(o), as.matrix(df_exclusive))
  expect_equal(colnames(o), colnames(df_exclusive))

  o2 <- with(df_exclusive, multi(a, b, c, inclusive = FALSE, type = "competing"))
  expect_s3_class(o2, "multi_competing")
  expect_equivalent(unclass(o2), as.matrix(df_exclusive))
  expect_equal(colnames(o2), colnames(df_exclusive))

  o3 <- with(df_exclusive, multi(A = a, B = b, C = c, inclusive = FALSE, type = "ordered"))
  expect_equal(colnames(o3), c("A", "B", "C"))

  o4 <- with(df_exclusive, multi(A = a, B = b, C = c, inclusive = FALSE, type = "competing"))
  expect_equal(colnames(o4), c("A", "B", "C"))

  # NA should be allowed in leftmost category for competing outcomes
  df_e2 <- dplyr::add_row(df_exclusive, a = NA, b = 5, c = 10)
  o5 <- with(df_e2, multi(a, b, c, inclusive = FALSE, type = "competing"))
  expect_equivalent(unclass(o5), as.matrix(df_e2))
})
