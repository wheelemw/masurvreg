#' @useDynLib rstanarm, .registration = TRUE 

.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  library(loo)
  library(ggplot2)
  for (m in modules) loadModule(m, what = TRUE)
}

.onAttach <- function(...) {
  rstanarmLib <- dirname(system.file(package = "rstanarm"))
  pkgdesc <- suppressWarnings(utils::packageDescription("rstanarm", lib.loc = rstanarmLib))
  if (length(pkgdesc) > 1) {
    builddate <- gsub(';.*$', '', pkgdesc$Packaged)
    packageStartupMessage(paste("survregstack (Version ", pkgdesc$Version, ", packaged: ", builddate, ")", sep = ""))
  }
  
  ggplot2::theme_set(bayesplot::theme_default())
}
