test_that("PostProcessChain() writes the grid files with the right shapes", {
  path <- gl_fixture_run()          # post-processed on a 20 x 20 grid
  ncell <- 20 * 20

  for (f in c("proba.pop.membership.txt", "proba.pop.membership.indiv.txt",
              "modal.pop.txt", "modal.pop.indiv.txt",
              "postprocess.parameters.txt", "perm.txt"))
    expect_true(file.exists(file.path(path, f)), info = f)

  dom <- as.matrix(utils::read.table(file.path(path, "proba.pop.membership.txt")))
  expect_equal(nrow(dom), ncell)
  expect_equal(ncol(dom), 2 + 5)    # x, y, then one column per npopmax

  modal <- as.matrix(utils::read.table(file.path(path, "modal.pop.txt")))
  expect_equal(nrow(modal), ncell)
  expect_equal(ncol(modal), 3)
})

test_that("membership probabilities are probabilities and sum to one", {
  path <- gl_fixture_run()
  dom <- as.matrix(utils::read.table(file.path(path, "proba.pop.membership.txt")))
  p <- dom[, -(1:2), drop = FALSE]

  expect_true(all(p >= 0 & p <= 1))
  expect_equal(unname(rowSums(p)), rep(1, nrow(p)), tolerance = 1e-6)
})

test_that("the modal cluster is the argmax of the membership probabilities", {
  path <- gl_fixture_run()
  dom <- as.matrix(utils::read.table(file.path(path, "proba.pop.membership.txt")))
  modal <- as.matrix(utils::read.table(file.path(path, "modal.pop.txt")))

  expect_equal(unname(modal[, 3]), unname(max.col(dom[, -(1:2), drop = FALSE])))
})

test_that("per-individual membership has one row per individual", {
  path <- gl_fixture_run()
  indiv <- as.matrix(utils::read.table(
    file.path(path, "proba.pop.membership.indiv.txt")))
  expect_equal(nrow(indiv), 100)
  p <- indiv[, -(1:2), drop = FALSE]
  expect_equal(unname(rowSums(p)), rep(1, 100), tolerance = 1e-6)
})

test_that("PostProcessChain() rejects a path that does not end in a separator", {
  dat <- gl_example_data()
  bad <- sub("/$", "", gl_fixture_run())
  expect_error(
    PostProcessChain(coordinates = dat$coord, path.mcmc = bad,
                     nxdom = 10, nydom = 10, burnin = 1))
})

test_that("a truncated chain file raises a catchable error, not a crash", {
  skip_on_cran()

  ## Before 5.0.0 the reads happened in Fortran with no iostat= guard, so a
  ## short file called abort() and killed the R process outright -- try()
  ## could not catch it. This is the regression guard for that.
  src <- gl_fixture_run_raw()
  dst <- gl_new_dir("truncated")
  file.copy(list.files(src, full.names = TRUE), dst)

  f <- file.path(dst, "frequencies.txt")
  writeLines(readLines(f, warn = FALSE)[1:5], f)

  dat <- gl_example_data()
  err <- tryCatch(
    utils::capture.output(
      PostProcessChain(coordinates = dat$coord, path.mcmc = dst,
                       nxdom = 10, nydom = 10, burnin = 1)),
    error = conditionMessage)

  expect_type(err, "character")
  expect_match(err, "frequencies.txt")
  expect_match(err, "truncated")
})

test_that("a missing chain file names the file and the write.* option", {
  skip_on_cran()

  src <- gl_fixture_run_raw()
  dst <- gl_new_dir("missing-chain")
  file.copy(list.files(src, full.names = TRUE), dst)
  file.remove(file.path(dst, "coord.nuclei.txt"))

  dat <- gl_example_data()
  err <- tryCatch(
    utils::capture.output(
      PostProcessChain(coordinates = dat$coord, path.mcmc = dst,
                       nxdom = 10, nydom = 10, burnin = 1)),
    error = conditionMessage)

  expect_match(err, "coord.nuclei.txt")
  expect_match(err, "missing")
})

test_that("no character argument is passed to .Fortran any more", {
  skip_on_cran()

  ## Passing character to .Fortran is deprecated and emits a warning; the
  ## whole point of moving the file reading into R was to stop doing it.
  ## Work on a copy: post-processing the shared raw fixture in place would
  ## break the tests that rely on it being un-post-processed.
  dat <- gl_example_data()
  path <- gl_new_dir("no-char-arg")
  file.copy(list.files(gl_fixture_run_raw(), full.names = TRUE), path)
  w <- character()
  withCallingHandlers(
    utils::capture.output(
      PostProcessChain(coordinates = dat$coord, path.mcmc = path,
                       nxdom = 10, nydom = 10, burnin = 1)),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    })
  expect_length(grep("char vector to .Fortran", w), 0)
})
