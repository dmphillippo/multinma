# Tests results of package vignettes

skip_on_cran()

test_that("Smoking example", {
  rmarkdown::render("../../vignettes/example_smoking.Rmd")
  file.remove("../../vignettes/example_smoking.html")
})

test_that("Thrombolytics example", {
  rmarkdown::render("../../vignettes/example_thrombolytics.Rmd")
  file.remove("../../vignettes/example_thrombolytics.html")
})

test_that("Blocker example", {
  rmarkdown::render("../../vignettes/example_blocker.Rmd")
  file.remove("../../vignettes/example_blocker.html")
})
