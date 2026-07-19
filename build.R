library(roxygen2)
#setwd("~/git/")
#devtools::create("analyticgrimu")
setwd("~/git/analyticgrimu")

devtools::document()

devtools::check(vignettes = FALSE)

#devtools::install(build_vignettes = TRUE)
#vignette("analyticgrimu")

# or from github, after push
devtools::install_github("ianhussey/analyticgrimu")

library(analyticgrimu)

?analyticgrimu

vignette("analyticgrimu")

detach("package:analyticgrimu", unload = TRUE)

# once you have the package updated, you can use it to build the vignettes, check the whole thing, and reinstall again
devtools::build_vignettes()
devtools::check()

# cran checks
# win-builder
library(devtools)
check_win_devel() # emails results to the maintainer address in DESCRIPTION

# R-hub
library(rhub)
# rhub_setup() # one time
rhub_check() # v2: runs on GitHub Actions in your repo
