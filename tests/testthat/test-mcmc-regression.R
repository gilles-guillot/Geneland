## Numerical regression on the sampler.
##
## These lock the sampler's output to known values so that refactoring the
## Fortran (in particular removing its file I/O, which is still pending)
## cannot silently change the numbers.
##
## They are skipped on CRAN on purpose. The chain's path depends on floating
## point accept/reject decisions, so a different compiler, optimisation level
## or BLAS can legitimately produce a different -- equally valid -- chain.
## Exact values are a regression guard for this codebase on a fixed toolchain,
## not a portable correctness claim. The portable invariants live in
## test-mcmc-validation.R and are not skipped.

test_that("the sampler is reproducible from a fixed seed", {
  skip_on_cran()

  dat <- gl_example_data()
  run_once <- function() {
    path <- gl_new_dir("repro")
    set.seed(4242)
    utils::capture.output(
      MCMC(coordinates = dat$coord, geno.dip.codom = dat$geno,
           path.mcmc = path, rate.max = nrow(dat$coord), delta.coord = 0,
           npopmin = 1, npopinit = 3, npopmax = 4,
           nb.nuclei.max = 3 * nrow(dat$coord),
           nit = 100, thinning = 10, freq.model = "Uncorrelated",
           varnpop = TRUE, spatial = TRUE, filter.null.alleles = FALSE))
    list(npop = scan(file.path(path, "populations.numbers.txt"), quiet = TRUE),
         ll   = scan(file.path(path, "log.likelihood.txt"), quiet = TRUE))
  }

  a <- run_once()
  b <- run_once()
  expect_identical(a$npop, b$npop)
  expect_identical(a$ll, b$ll)
})

test_that("sampler output matches the recorded reference values", {
  skip_on_cran()

  ## Reference produced by Geneland 5.0.0 on x86_64-pc-linux-gnu, gfortran 11.
  ## These are byte-identical to the values produced by version 4.9.2, which
  ## is the evidence that the 5.0.0 refactoring left the sampler alone.
  expected_npop <- c(3, 3, 3, rep(2, 27))
  expected_npp <- c(52, 52, 51, 52, 51, 49, 45, 45, 47, 46, 45, 43, 44, 43,
                    44, 42, 42, 44, 43, 45, 46, 43, 45, 47, 40, 45, 43, 41,
                    42, 42)
  expected_ll_head <- c(-6098.9555147403, -6017.0915527358, -5941.5102980563)

  path <- gl_fixture_run()
  npop <- scan(file.path(path, "populations.numbers.txt"), quiet = TRUE)
  npp  <- scan(file.path(path, "nuclei.numbers.txt"), quiet = TRUE)
  ll   <- scan(file.path(path, "log.likelihood.txt"), quiet = TRUE)

  expect_equal(npop, expected_npop)
  expect_equal(npp, expected_npp)
  expect_equal(ll[1:3], expected_ll_head, tolerance = 1e-8)
})

test_that("the Correlated frequency model also runs and writes drifts", {
  skip_on_cran()

  dat <- gl_example_data()
  path <- gl_new_dir("correlated")
  set.seed(11)
  utils::capture.output(
    MCMC(coordinates = dat$coord, geno.dip.codom = dat$geno, path.mcmc = path,
         rate.max = nrow(dat$coord), delta.coord = 0,
         npopmin = 1, npopinit = 3, npopmax = 4,
         nb.nuclei.max = 3 * nrow(dat$coord),
         nit = 100, thinning = 10, freq.model = "Correlated",
         varnpop = TRUE, spatial = TRUE, filter.null.alleles = FALSE))

  drift <- as.matrix(utils::read.table(file.path(path, "drifts.txt")))
  expect_equal(dim(drift), c(10L, 4L))

  ## Clusters absent from the model at a given iteration are written as the
  ## sentinel -999; every real drift factor is a probability.
  real <- drift[drift != -999]
  expect_true(length(real) > 0)
  expect_true(all(real >= 0 & real <= 1))
  expect_true(all(drift == -999 | (drift >= 0 & drift <= 1)))
})

test_that("autoplot('drift') drops the -999 sentinel instead of plotting it", {
  skip_on_cran()

  dat <- gl_example_data()
  path <- gl_new_dir("drift-sentinel")
  set.seed(11)
  utils::capture.output(
    MCMC(coordinates = dat$coord, geno.dip.codom = dat$geno, path.mcmc = path,
         rate.max = nrow(dat$coord), delta.coord = 0,
         npopmin = 1, npopinit = 3, npopmax = 4,
         nb.nuclei.max = 3 * nrow(dat$coord),
         nit = 100, thinning = 10, freq.model = "Correlated",
         varnpop = TRUE, spatial = TRUE, filter.null.alleles = FALSE))

  raw <- as.matrix(utils::read.table(file.path(path, "drifts.txt")))
  expect_true(any(raw == -999))          # the fixture really does contain it

  d <- ggplot2::ggplot_build(autoplot(read_geneland(path), "drift",
                                      burnin = 0))$data[[1]]
  plotted <- d$y[!is.na(d$y)]
  expect_true(all(plotted >= 0 & plotted <= 1))
  expect_false(any(plotted < 0))
})
