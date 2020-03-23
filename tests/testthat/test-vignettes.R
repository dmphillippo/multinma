# Tests results of package vignettes

skip_on_cran()

test_that("Smoking example", {
  rmarkdown::render("../../vignettes/example_smoking.Rmd")
  file.remove("../../vignettes/example_smoking.html")
})
