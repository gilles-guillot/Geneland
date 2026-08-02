test_that("simdata(IBD = TRUE) explains how to get RandomFields when it is absent", {
  skip_if(requireNamespace("RandomFields", quietly = TRUE),
          "RandomFields is installed, so the guard cannot be exercised")

  ## RandomFields was archived from CRAN, so this path must fail with
  ## instructions rather than 'could not find function "GaussRF"'.
  err <- tryCatch(
    simdata(nindiv = 20, coord.lim = c(0, 1, 0, 1), npop = 2,
            number.nuclei = 10, allele.numbers = rep(3, 4),
            sim.gen = TRUE, IBD = TRUE, model = "whittle",
            beta = 0.1, gamma = 1, seed.freq = 1),
    error = conditionMessage)

  expect_match(err, "RandomFields")
  expect_match(err, "install_github")
  expect_match(err, "IBD = FALSE")
})

test_that("show.simdata() asks for fields only when it needs it", {
  skip_if(requireNamespace("fields", quietly = TRUE),
          "fields is installed, so the guard cannot be exercised")
  skip_on_cran()

  set.seed(1)
  d <- simdata(nindiv = 20, coord.lim = c(0, 1, 0, 1), npop = 2,
               number.nuclei = 10, allele.numbers = rep(3, 4),
               sim.gen = TRUE, IBD = FALSE)
  expect_error(
    show.simdata(d, plot.freq.indiv = TRUE, file.plot.freq.indiv = NA),
    "fields")
})

test_that("run_geneland_app() reports which interface packages are missing", {
  skip_if(all(vapply(c("shiny", "bslib"), requireNamespace,
                     logical(1), quietly = TRUE)),
          "shiny and bslib are installed")
  expect_error(run_geneland_app(), "install.packages")
})

test_that("the Shiny application is installed with the package", {
  app <- system.file("shiny", "app.R", package = "Geneland")
  expect_true(nzchar(app))
  expect_true(file.exists(app))
  expect_no_error(parse(app))
})
