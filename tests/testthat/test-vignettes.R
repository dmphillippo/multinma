# Tests results of package vignettes

skip_on_cran()

td <- tempdir()

test_that("Smoking example", {
  rmarkdown::render("../../vignettes/example_smoking.Rmd", output_dir = td)
})

test_that("Thrombolytics example", {
  rmarkdown::render("../../vignettes/example_thrombolytics.Rmd", output_dir = td)
})

test_that("Blocker example", {
  rmarkdown::render("../../vignettes/example_blocker.Rmd", output_dir = td)
})

test_that("Dietary fat example", {
  rmarkdown::render("../../vignettes/example_dietary_fat.Rmd", output_dir = td)
})

test_that("Diabetes example", {
  rmarkdown::render("../../vignettes/example_diabetes.Rmd", output_dir = td)
})

test_that("Parkinsons example", {
  rmarkdown::render("../../vignettes/example_parkinsons.Rmd", output_dir = td)
})


