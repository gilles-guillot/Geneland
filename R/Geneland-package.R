#' Geneland: clustering of georeferenced genetic and morphometric data
#' 
#' Geneland implements Bayesian models and Markov chain Monte Carlo (MCMC)
#' algorithms to infer population genetic structure from georeferenced
#' individual multi-locus genotypes and quantitative phenotypes.
#' 
#' @section Getting started:
#' [MCMC()] runs the sampler and writes its output to a directory,
#' [PostProcessChain()] post-processes a chain to deal with label switching,
#' and the `Plot*` functions and [PosteriorMode()] display the results. See
#' `vignette("Geneland")` for a worked example.
#' 
#' @section Models for allele frequencies:
#' The `freq.model` argument of [MCMC()] selects how allele frequencies are
#' treated across populations.
#' 
#' \describe{
#'   \item{`"Uncorrelated"`}{Frequencies in each population are given
#'     independent Dirichlet priors. This is the safer default: it does not
#'     assume the populations share a common history, and it is the only
#'     model compatible with `filter.null.alleles = TRUE`.}
#'   \item{`"Correlated"`}{Frequencies are drawn around those of a common
#'     ancestral population, each population having its own drift parameter
#'     with a `Beta(shape1, shape2)` prior. This model has more power to
#'     detect subtle structure, but it is also more prone to detecting
#'     spurious clusters when its assumptions are violated.}
#' }
#' 
#' @keywords internal
#' @aliases Geneland-package Geneland
"_PACKAGE"

## Namespace imports -----------------------------------------------------
## Base-package functions used throughout the package. Declared in one place
## rather than scattered across files so the full set stays visible.

#' @importFrom graphics abline contour hist image lines par plot.default points polygon text title
#' @importFrom grDevices dev.new dev.off heat.colors palette pdf postscript terrain.colors
#' @importFrom parallel clusterApply makeCluster stopCluster
#' @importFrom rlang .data
#' @importFrom stats dist lm pnorm qexp qgamma rbeta rexp rgamma rnorm rpois runif
#' @importFrom utils browseURL file_test help.start installed.packages packageVersion read.table write.table
##
## Geneland.GUI uses a large number of tcltk functions; import the whole
## namespace rather than maintain a list that silently rots.
#' @import tcltk
##
#' @useDynLib Geneland
NULL
