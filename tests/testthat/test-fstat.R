test_that("Fstat() returns Fis and a pairwise Fst matrix", {
  dat <- gl_example_data()
  fs <- Fstat(genotypes = dat$geno, npop = 2,
              pop.mbrship = rep(1:2, each = 50))

  expect_named(fs, c("Fis", "Fst"))
  expect_length(fs$Fis, 2)
  expect_equal(dim(fs$Fst), c(2L, 2L))
})

test_that("the pairwise Fst matrix is symmetric with a zero diagonal", {
  dat <- gl_example_data()
  fs <- Fstat(genotypes = dat$geno, npop = 3,
              pop.mbrship = rep(1:3, length.out = 100))

  expect_equal(diag(fs$Fst), rep(0, 3))
  lower <- fs$Fst[lower.tri(fs$Fst)]
  upper <- t(fs$Fst)[lower.tri(fs$Fst)]
  keep <- !is.na(lower) & !is.na(upper)
  expect_equal(lower[keep], upper[keep])
})

test_that("F statistics fall in their admissible range", {
  dat <- gl_example_data()
  fs <- Fstat(genotypes = dat$geno, npop = 2,
              pop.mbrship = rep(1:2, each = 50))

  fis <- fs$Fis[!is.na(fs$Fis)]
  fst <- fs$Fst[!is.na(fs$Fst)]
  expect_true(all(fis >= -1 & fis <= 1))
  expect_true(all(fst >= -1 & fst <= 1))
})

test_that("a random partition of a panmictic sample gives near-zero Fst", {
  ## Splitting one population at random creates no real differentiation, so
  ## the estimator should sit near zero. A loose bound: this is a property
  ## check, not a claim about the estimator's variance.
  dat <- gl_example_data()
  set.seed(1)
  fs <- Fstat(genotypes = dat$geno, npop = 2,
              pop.mbrship = sample(rep(1:2, each = 50)))
  expect_lt(abs(fs$Fst[1, 2]), 0.05)
})

test_that("Fstat() matches its recorded values", {
  skip_on_cran()
  dat <- gl_example_data()
  fs <- Fstat(genotypes = dat$geno, npop = 2,
              pop.mbrship = rep(1:2, each = 50))
  ## testthat's tolerance is relative, so small quantities need full digits.
  expect_equal(fs$Fis, c(0.046139091195, 0.055886400776), tolerance = 1e-6)
  expect_equal(fs$Fst[1, 2], 0.000390022446295, tolerance = 1e-6)
})

test_that("Fstat() refuses haploid data rather than silently misbehaving", {
  dat <- gl_example_data()
  expect_error(
    Fstat(genotypes = dat$geno, npop = 2,
          pop.mbrship = rep(1:2, each = 50), ploidy = 1),
    "haploid")
})
