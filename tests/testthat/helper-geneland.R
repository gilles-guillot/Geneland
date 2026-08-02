## Shared fixtures.
##
## MCMC() is the slowest thing in the suite, so a single short chain is run
## once per test session and reused. Everything is written under tempdir().

gl_example_data <- function() {
  d <- system.file("extdata", package = "Geneland")
  list(
    coord = as.matrix(utils::read.table(file.path(d, "coordinates.txt"))),
    geno  = as.matrix(utils::read.table(file.path(d, "genotypes.txt")))
  )
}

## A path ending in the separator, as MCMC() requires.
gl_new_dir <- function(name) {
  p <- file.path(tempdir(), paste0(name, "-", as.integer(runif(1, 1e6, 9e6))))
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  paste0(p, "/")
}

.gl_cache <- new.env(parent = emptyenv())

## A short chain plus its post-processing. Deliberately tiny: the point is to
## exercise the code paths, not to sample the posterior properly.
gl_fixture_run <- function() {
  if (!is.null(.gl_cache$path)) return(.gl_cache$path)

  dat <- gl_example_data()
  path <- gl_new_dir("fixture")
  set.seed(20240101)
  utils::capture.output(
    MCMC(coordinates = dat$coord, geno.dip.codom = dat$geno, path.mcmc = path,
         rate.max = nrow(dat$coord), delta.coord = 0,
         npopmin = 1, npopinit = 3, npopmax = 5,
         nb.nuclei.max = 3 * nrow(dat$coord),
         nit = 300, thinning = 10, freq.model = "Uncorrelated",
         varnpop = TRUE, spatial = TRUE, filter.null.alleles = FALSE))
  utils::capture.output(
    PostProcessChain(coordinates = dat$coord, path.mcmc = path,
                     nxdom = 20, nydom = 20, burnin = 5))

  .gl_cache$path <- path
  path
}

## A chain only, without post-processing, for testing the "you must
## post-process first" paths.
gl_fixture_run_raw <- function() {
  if (!is.null(.gl_cache$raw)) return(.gl_cache$raw)

  dat <- gl_example_data()
  path <- gl_new_dir("fixture-raw")
  set.seed(99)
  utils::capture.output(
    MCMC(coordinates = dat$coord, geno.dip.codom = dat$geno, path.mcmc = path,
         rate.max = nrow(dat$coord), delta.coord = 0,
         npopmin = 1, npopinit = 2, npopmax = 3,
         nb.nuclei.max = 3 * nrow(dat$coord),
         nit = 100, thinning = 10, freq.model = "Uncorrelated",
         varnpop = TRUE, spatial = TRUE, filter.null.alleles = FALSE))

  .gl_cache$raw <- path
  path
}
