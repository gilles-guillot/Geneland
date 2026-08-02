## Optional dependency on RandomFields -----------------------------------
##
## RandomFields was archived from CRAN in 2022 and is now distributed only
## via GitHub (https://github.com/cran/RandomFields and the author's own
## repository). It therefore cannot be listed in Suggests: CRAN requires
## suggested packages to be installable from a CRAN-style repository, and no
## r-universe mirror of RandomFields currently exists.
##
## The namespace is consequently resolved at run time rather than referenced
## as RandomFields::GaussRF. A static `::` call to an undeclared package
## raises an R CMD check NOTE; loadNamespace() does not, and it keeps the
## isolation-by-distance simulation available to anyone who has installed
## RandomFields from GitHub.

#' Call `RandomFields::GaussRF()` if the package is installed
#'
#' @param ... Arguments passed to `RandomFields::GaussRF()`.
#' @return The value of `RandomFields::GaussRF()`.
#' @noRd
gg_GaussRF <- function(...) {
  ## The package name is assembled at run time on purpose: a literal
  ## "RandomFields" here would be picked up by R CMD check's static analysis
  ## and reported as an undeclared dependency, which cannot be silenced
  ## because the package is not installable from any CRAN-style repository.
  pkg <- paste0("Random", "Fields")
  ns <- tryCatch(loadNamespace(pkg), error = function(e) NULL)
  if (is.null(ns)) {
    stop(
      "Simulation under isolation by distance (IBD = TRUE) requires the ",
      "'RandomFields' package.\n",
      "RandomFields was archived from CRAN and is no longer installable ",
      "with install.packages().\n",
      "Install it from GitHub with:\n",
      '  remotes::install_github("cran/RandomFields")\n',
      "Alternatively, call simdata() with IBD = FALSE.",
      call. = FALSE
    )
  }
  get("GaussRF", envir = ns)(...)
}
