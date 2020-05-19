# Precompile vignettes that run Stan models

root <- rprojroot::is_r_package

rmarkdown::render(root$find_file("vignettes", "example_smoking.Rmd"))
rmarkdown::render(root$find_file("vignettes", "example_thrombolytics.Rmd"))
rmarkdown::render(root$find_file("vignettes", "example_blocker.Rmd"))
rmarkdown::render(root$find_file("vignettes", "example_dietary_fat.Rmd"))
rmarkdown::render(root$find_file("vignettes", "example_diabetes.Rmd"))
rmarkdown::render(root$find_file("vignettes", "example_parkinsons.Rmd"))
rmarkdown::render(root$find_file("vignettes", "example_transfusion.Rmd"))
rmarkdown::render(root$find_file("vignettes", "example_atrial_fibrillation.Rmd"))
rmarkdown::render(root$find_file("vignettes", "example_statins.Rmd"))
rmarkdown::render(root$find_file("vignettes", "example_bcg_vaccine.Rmd"))
rmarkdown::render(root$find_file("vignettes", "example_plaque_psoriasis.Rmd"))
