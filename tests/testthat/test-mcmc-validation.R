## Portable properties of MCMC(): argument validation and invariants that hold
## whatever the toolchain. These are NOT skipped on CRAN.

test_that("MCMC() tolerates a path without a trailing separator", {
  ## The documentation says the path "has to" end with a slash, but the code
  ## appends one. Pin the forgiving behaviour so it is not lost by accident.
  skip_on_cran()
  dat <- gl_example_data()
  bare <- file.path(tempdir(), paste0("no-slash-", as.integer(runif(1, 1e6, 9e6))))
  dir.create(bare, showWarnings = FALSE)

  utils::capture.output(
    MCMC(coordinates = dat$coord, geno.dip.codom = dat$geno, path.mcmc = bare,
         rate.max = nrow(dat$coord), npopmin = 1, npopinit = 2, npopmax = 3,
         nb.nuclei.max = 300, nit = 20, thinning = 10,
         freq.model = "Uncorrelated", filter.null.alleles = FALSE))

  ## Output must land inside the directory, not beside it with a glued name.
  expect_true(file.exists(file.path(bare, "parameters.txt")))
  expect_false(file.exists(paste0(bare, "parameters.txt")))
})

test_that("MCMC() rejects a directory that does not exist", {
  dat <- gl_example_data()
  missing <- paste0(file.path(tempdir(), "definitely-not-here-12345"), "/")
  expect_error(
    MCMC(coordinates = dat$coord, geno.dip.codom = dat$geno,
         path.mcmc = missing,
         rate.max = 100, npopmin = 1, npopinit = 2, npopmax = 3,
         nb.nuclei.max = 300, nit = 10, thinning = 1),
    regexp = "[Dd]irectory")
})

test_that("dominant diploid and haploid data cannot be combined", {
  dat <- gl_example_data()
  dom <- matrix(0L, nrow = nrow(dat$geno), ncol = 5)
  hap <- matrix(1L, nrow = nrow(dat$geno), ncol = 5)
  expect_error(
    MCMC(coordinates = dat$coord, geno.dip.dom = dom, geno.hap = hap,
         path.mcmc = gl_new_dir("combo"),
         rate.max = 100, npopmin = 1, npopinit = 2, npopmax = 3,
         nb.nuclei.max = 300, nit = 10, thinning = 1),
    regexp = "jointly")
})

test_that("a chain writes the expected files with the expected shapes", {
  path <- gl_fixture_run()
  nsaved <- 300 / 10

  for (f in c("parameters.txt", "populations.numbers.txt",
              "nuclei.numbers.txt", "log.likelihood.txt",
              "log.posterior.density.txt", "coord.nuclei.txt",
              "color.nuclei.txt"))
    expect_true(file.exists(file.path(path, f)), info = f)

  npop <- scan(file.path(path, "populations.numbers.txt"), quiet = TRUE)
  npp  <- scan(file.path(path, "nuclei.numbers.txt"), quiet = TRUE)
  ll   <- scan(file.path(path, "log.likelihood.txt"), quiet = TRUE)

  expect_length(npop, nsaved)
  expect_length(npp, nsaved)
  expect_length(ll, nsaved)
})

test_that("the sampled quantities stay inside the bounds that were asked for", {
  path <- gl_fixture_run()
  npop <- scan(file.path(path, "populations.numbers.txt"), quiet = TRUE)
  npp  <- scan(file.path(path, "nuclei.numbers.txt"), quiet = TRUE)

  ## npopmin = 1, npopmax = 5, nb.nuclei.max = 300 in the fixture.
  expect_true(all(npop >= 1), info = "number of clusters below npopmin")
  expect_true(all(npop <= 5), info = "number of clusters above npopmax")
  expect_true(all(npop == round(npop)), info = "cluster counts must be whole")
  expect_true(all(npp >= 1 & npp <= 300))

  ll <- scan(file.path(path, "log.likelihood.txt"), quiet = TRUE)
  expect_true(all(is.finite(ll)), info = "non-finite log-likelihood")
  expect_true(all(ll < 0), info = "log-likelihood of discrete data must be < 0")
})

test_that("parameters.txt round-trips the settings that were used", {
  path <- gl_fixture_run()
  p <- utils::read.table(file.path(path, "parameters.txt"),
                         stringsAsFactors = FALSE)
  get <- function(k) p[[3]][p[[1]] == k]

  expect_equal(as.numeric(get("nit")), 300)
  expect_equal(as.numeric(get("thinning")), 10)
  expect_equal(as.numeric(get("npopmax")), 5)
  expect_equal(as.numeric(get("nindiv")), 100)
  expect_equal(get("freq.model"), "Uncorrelated")
})
