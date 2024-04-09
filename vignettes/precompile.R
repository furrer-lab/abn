# Pre-compile vignettes that depend on the INLA package
knitr::knit("vignettes/data_simulation.Rmd.orig", "vignettes/data_simulation.Rmd")
knitr::knit("vignettes/parameter_learning.Rmd.orig", "vignettes/parameter_learning.Rmd")
knitr::knit("vignettes/multiprocessing.Rmd.orig", "vignettes/multiprocessing.Rmd")
