library(dplyr)
library(tidyr)

test_that("prefer_n argument", {
  m <- "vector of integers greater than 1"
  expect_error(breaks_integer(prefer_n = "a"), m)
  expect_error(breaks_integer(prefer_n = -1), m)
  expect_error(breaks_integer(prefer_n = 0), m)
  expect_error(breaks_integer(prefer_n = 1:2), m)
  expect_error(breaks_integer(prefer_n = 1.5), m)
  expect_error(breaks_integer(prefer_n = Inf), m)
  expect_error(breaks_integer(prefer_n = list()), m)
  expect_error(breaks_integer(prefer_n = NA), m)
  expect_error(breaks_integer(prefer_n = NULL), m)
})

test_that("extend argument", {
  m <- "logical value TRUE or FALSE"
  expect_error(breaks_integer(extend = "a"), m)
  expect_error(breaks_integer(extend = 1), m)
  expect_error(breaks_integer(extend = list()), m)
  expect_error(breaks_integer(extend = NA), m)
  expect_error(breaks_integer(extend = NULL), m)
  expect_error(breaks_integer(extend = c(TRUE, FALSE)), m)
})

test_that("positive argument", {
  m <- "logical value TRUE or FALSE"
  expect_error(breaks_integer(positive = "a"), m)
  expect_error(breaks_integer(positive = 1), m)
  expect_error(breaks_integer(positive = list()), m)
  expect_error(breaks_integer(positive = NA), m)
  expect_error(breaks_integer(positive = NULL), m)
  expect_error(breaks_integer(positive = c(TRUE, FALSE)), m)
})

bdat1 <- tibble(x0 = 1:10001) %>%
  rowwise() %>%
  mutate(x = list(1:x0),
         b_int = list(breaks_integer(extend = FALSE, positive = TRUE)(x)),
         b_int_ext = list(breaks_integer(extend = TRUE, positive = TRUE)(x))
  )

bdat2 <- tibble(x0 = 2:10001) %>%
  rowwise() %>%
  mutate(x = list(2:x0),
         b_int = list(breaks_integer(extend = FALSE, positive = TRUE)(x)),
         b_int_ext = list(breaks_integer(extend = TRUE, positive = TRUE)(x))
)

bdat5 <- tibble(x0 = 5:10001) %>%
  rowwise() %>%
  mutate(x = list(5:x0),
         b_int = list(breaks_integer(extend = FALSE, positive = TRUE)(x)),
         b_int_ext = list(breaks_integer(extend = TRUE, positive = TRUE)(x))
)

bdat_neg <- tibble(x0 = -100:100) %>%
  rowwise() %>%
  mutate(x = list(-100:x0),
         b_int = list(breaks_integer(extend = FALSE, positive = FALSE)(x)),
         b_int_ext = list(breaks_integer(extend = TRUE, positive = FALSE)(x))
  )

test_that("Breaks are equally spaced with extend = TRUE", {
  expect_true(all(purrr::map_lgl(bdat1$b_int_ext[-1], ~n_distinct(diff(.)) == 1)))
  expect_true(all(purrr::map_lgl(bdat2$b_int_ext[-1], ~n_distinct(diff(.)) == 1)))
  expect_true(all(purrr::map_lgl(bdat5$b_int_ext[-1], ~n_distinct(diff(.)) == 1)))
  expect_true(all(purrr::map_lgl(bdat_neg$b_int_ext[-1], ~n_distinct(diff(.)) == 1)))
})

test_that("Breaks are integers, and > 0 when positive = TRUE", {
  expect_true(all(purrr::map_lgl(bdat1$b_int_ext, ~all(rlang::is_integerish(.)) && all(. > 0))))
  expect_true(all(purrr::map_lgl(bdat2$b_int_ext, ~all(rlang::is_integerish(.)) && all(. > 0))))
  expect_true(all(purrr::map_lgl(bdat5$b_int_ext, ~all(rlang::is_integerish(.)) && all(. > 0))))
  expect_true(all(purrr::map_lgl(bdat1$b_int, ~all(rlang::is_integerish(.)) && all(. > 0))))
  expect_true(all(purrr::map_lgl(bdat2$b_int, ~all(rlang::is_integerish(.)) && all(. > 0))))
  expect_true(all(purrr::map_lgl(bdat5$b_int, ~all(rlang::is_integerish(.)) && all(. > 0))))
  expect_true(all(purrr::map_lgl(bdat_neg$b_int, ~all(rlang::is_integerish(.)))))
  expect_true(all(purrr::map_lgl(bdat_neg$b_int_ext, ~all(rlang::is_integerish(.)))))
})
