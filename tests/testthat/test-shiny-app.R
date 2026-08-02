## The Shiny interface, driven headlessly through shiny::testServer.
## Skipped on CRAN: it runs a (short) chain and needs the optional interface
## packages.

test_that("the app runs the whole workflow: load, sample, post-process, plot", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")

  app <- shiny::shinyAppFile(
    system.file("shiny", "app.R", package = "Geneland"))

  shiny::testServer(app, {
    session$setInputs(src = "example", load = 1)
    expect_equal(nrow(rv$coord), 100)
    expect_equal(ncol(rv$geno) / 2, 20)

    session$setInputs(freq_model = "Uncorrelated", spatial = TRUE,
                      varnpop = TRUE, filterna = FALSE,
                      npopmin = 1, npopinit = 3, npopmax = 4,
                      nit = 200, thinning = 10, delta_coord = 0)
    expect_equal(derived()$nsaved, 20)
    expect_equal(derived()$rate_max, 100)
    expect_equal(derived()$nb_nuclei_max, 300)

    session$setInputs(run = 1)
    expect_false(is.null(rv$run))
    expect_s3_class(rv$run, "geneland_run")
    expect_match(rv$status, "MCMC finished")

    session$setInputs(burnin = 5, nxdom = 15, nydom = 15, post = 1)
    expect_true(rv$posted)

    for (ty in c("npop", "npop_post", "loglik", "rate")) {
      session$setInputs(diag = ty, diag_burnin = 5)
      expect_s3_class(diag_plot(), "ggplot")
    }
    for (ty in c("map", "proba")) {
      session$setInputs(maptype = ty, show_pts = TRUE)
      expect_s3_class(map_plot(), "ggplot")
    }
  })
})

test_that("the app refuses mismatched inputs instead of failing later", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")

  app <- shiny::shinyAppFile(
    system.file("shiny", "app.R", package = "Geneland"))

  shiny::testServer(app, {
    dat <- gl_example_data()
    d <- file.path(tempdir(), "shiny-bad"); dir.create(d, showWarnings = FALSE)

    ## Coordinates with the wrong number of rows.
    fc <- file.path(d, "coord.txt"); fg <- file.path(d, "geno.txt")
    utils::write.table(dat$coord[1:50, ], fc, row.names = FALSE, col.names = FALSE)
    utils::write.table(dat$geno, fg, row.names = FALSE, col.names = FALSE)

    session$setInputs(
      src = "upload", header = FALSE,
      f_coord = list(datapath = fc, name = "coord.txt"),
      f_geno  = list(datapath = fg, name = "geno.txt"),
      load = 1)

    ## The bad upload must not become the active dataset.
    expect_null(rv$coord)
  })
})

test_that("the app leaves the user's filespace alone", {
  skip_on_cran()
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")

  app <- shiny::shinyAppFile(
    system.file("shiny", "app.R", package = "Geneland"))

  shiny::testServer(app, {
    session$setInputs(src = "example", load = 1)
    session$setInputs(freq_model = "Uncorrelated", spatial = TRUE,
                      varnpop = TRUE, filterna = FALSE,
                      npopmin = 1, npopinit = 2, npopmax = 3,
                      nit = 100, thinning = 10, delta_coord = 0)
    session$setInputs(run = 1)
    ## Everything must land under tempdir().
    expect_true(startsWith(normalizePath(rv$path),
                           normalizePath(tempdir())))
  })
})
