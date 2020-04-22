# Tests results of package vignettes

skip_on_cran()

td <- tempdir()

test_that("Smoking example", {
  rmarkdown::render("../../vignettes/example_smoking.Rmd", output_dir = td, params = list(run_tests = TRUE))
})

test_that("Thrombolytics example", {
  rmarkdown::render("../../vignettes/example_thrombolytics.Rmd", output_dir = td, params = list(run_tests = TRUE))
})

test_that("Blocker example", {
  rmarkdown::render("../../vignettes/example_blocker.Rmd", output_dir = td, params = list(run_tests = TRUE))
})

test_that("Dietary fat example", {
  rmarkdown::render("../../vignettes/example_dietary_fat.Rmd", output_dir = td, params = list(run_tests = TRUE))
})

test_that("Diabetes example", {
  rmarkdown::render("../../vignettes/example_diabetes.Rmd", output_dir = td, params = list(run_tests = TRUE))
})

test_that("Parkinsons example", {
  rmarkdown::render("../../vignettes/example_parkinsons.Rmd", output_dir = td, params = list(run_tests = TRUE))
})

test_that("Transfusion example", {
  rmarkdown::render("../../vignettes/example_transfusion.Rmd", output_dir = td, params = list(run_tests = TRUE))
})

test_that("Atrial fibrillation example", {
  rmarkdown::render("../../vignettes/example_atrial_fibrillation.Rmd", output_dir = td, params = list(run_tests = TRUE))
})

test_that("Statins example", {
  rmarkdown::render("../../vignettes/example_statins.Rmd", output_dir = td, params = list(run_tests = TRUE))
})

test_that("BCG vaccine example", {
  rmarkdown::render("../../vignettes/example_bcg_vaccine.Rmd", output_dir = td, params = list(run_tests = TRUE))
})

