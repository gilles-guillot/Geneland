test_that("every chain diagnostic builds", {
  run <- read_geneland(gl_fixture_run())
  for (ty in c("npop", "npop_post", "ntile", "ntile_post",
               "loglik", "logpost", "rate")) {
    p <- autoplot(run, ty, burnin = 5)
    expect_s3_class(p, "ggplot")
    ## Building is what catches bad aesthetics; constructing the object alone
    ## does not evaluate them.
    expect_no_error(ggplot2::ggplot_build(p))
  }
})

test_that("the maps build, with and without individual coordinates", {
  run <- read_geneland(gl_fixture_run())
  dat <- gl_example_data()
  for (ty in c("map", "proba")) {
    expect_no_error(ggplot2::ggplot_build(autoplot(run, ty)))
    expect_no_error(
      ggplot2::ggplot_build(autoplot(run, ty, coordinates = dat$coord)))
  }
})

test_that("maps require post-processing and say so", {
  run <- read_geneland(gl_fixture_run_raw())
  expect_error(autoplot(run, "map"), "PostProcessChain")
  expect_error(autoplot(run, "proba"), "PostProcessChain")
})

test_that("an unknown display type is rejected", {
  run <- read_geneland(gl_fixture_run())
  expect_error(autoplot(run, "not-a-type"))
})

test_that("burnin is validated", {
  run <- read_geneland(gl_fixture_run())
  expect_error(autoplot(run, "npop", burnin = -1), "non-negative")
  expect_error(autoplot(run, "npop", burnin = "ten"), "non-negative")
  ## The fixture has 30 saved iterations.
  expect_error(autoplot(run, "npop", burnin = 30), "smaller than")
  expect_error(autoplot(run, "npop", burnin = 999), "smaller than")
  expect_s3_class(autoplot(run, "npop", burnin = 0), "ggplot")
})

test_that("burnin defaults to the value used by PostProcessChain", {
  run <- read_geneland(gl_fixture_run())
  expect_equal(run$postprocess$burnin, 5)
  expect_no_error(ggplot2::ggplot_build(autoplot(run, "npop")))
})

test_that("coordinates must have two columns", {
  run <- read_geneland(gl_fixture_run())
  expect_error(autoplot(run, "map", coordinates = matrix(1:9, ncol = 3)),
               "two columns")
})

test_that("the posterior of the number of clusters spans the whole prior support", {
  run <- read_geneland(gl_fixture_run())
  d <- ggplot2::ggplot_build(autoplot(run, "npop_post", burnin = 5))$data[[1]]
  ## npopmin = 1, npopmax = 5 -> one bar per supported value, so a degenerate
  ## posterior cannot render as a single bar filling the panel.
  expect_equal(sort(unique(d$x)), 1:5)
  expect_equal(sum(d$y), 1, tolerance = 1e-8)
})

test_that("the membership map labels every cluster it draws", {
  run <- read_geneland(gl_fixture_run())
  b <- ggplot2::ggplot_build(autoplot(run, "map"))
  raster <- b$data[[1]]
  labels <- b$data[[length(b$data)]]
  ## Colour alone cannot carry identity on a map, so each region is labelled.
  expect_equal(nrow(labels), length(unique(raster$fill)))
})

test_that("probability surfaces are faceted, one panel per cluster, on [0,1]", {
  run <- read_geneland(gl_fixture_run())
  p <- autoplot(run, "proba")
  b <- ggplot2::ggplot_build(p)
  npop <- read_geneland(gl_fixture_run())$params$npopmax
  expect_equal(length(unique(b$data[[1]]$PANEL)), npop)
  expect_equal(p$scales$get_scales("fill")$limits, c(0, 1))
})

test_that("plot() prints and returns its input invisibly", {
  run <- read_geneland(gl_fixture_run())
  pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(run, type = "npop", burnin = 5))
})

test_that("plotting does not disturb the caller's graphics state", {
  run <- read_geneland(gl_fixture_run())
  pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  before <- graphics::par(c("mfrow", "mar"))
  ndev <- length(grDevices::dev.list())
  print(autoplot(run, "npop", burnin = 5))
  expect_identical(graphics::par(c("mfrow", "mar")), before)
  ## No stray dev.new(), unlike the base-graphics functions.
  expect_equal(length(grDevices::dev.list()), ndev)
})
