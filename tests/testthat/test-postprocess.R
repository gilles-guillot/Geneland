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
