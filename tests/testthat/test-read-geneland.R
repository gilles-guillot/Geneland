test_that("read_geneland() rejects bad input", {
  expect_error(read_geneland(123), "single character string")
  expect_error(read_geneland(c("a", "b")), "single character string")
  expect_error(read_geneland(file.path(tempdir(), "nope-98765")),
               "not found")

  empty <- file.path(tempdir(), "empty-dir-for-test")
  dir.create(empty, showWarnings = FALSE)
  expect_error(read_geneland(empty), "parameters.txt")
})

test_that("read_geneland() parses the run parameters", {
  run <- read_geneland(gl_fixture_run())

  expect_s3_class(run, "geneland_run")
  expect_equal(run$params$nit, 300)
  expect_equal(run$params$thinning, 10)
  expect_equal(run$params$npopmax, 5)
  expect_equal(run$params$nindiv, 100)

  ## Types are coerced, not left as character.
  expect_type(run$params$nit, "double")
  expect_type(run$params$spatial, "logical")
  expect_true(run$params$spatial)
  expect_type(run$params$freq.model, "character")
})

test_that("post-processing parameters are read when present and NULL when not", {
  expect_null(read_geneland(gl_fixture_run_raw())$postprocess)

  pp <- read_geneland(gl_fixture_run())$postprocess
  expect_false(is.null(pp))
  expect_equal(pp$nxdom, 20)
  expect_equal(pp$nydom, 20)
})

test_that("print() reports the run without error", {
  out <- utils::capture.output(print(read_geneland(gl_fixture_run())))
  expect_true(any(grepl("geneland_run", out)))
  expect_true(any(grepl("individuals", out)))
  expect_true(any(grepl("200 saved|30 saved", out)))
})

test_that("print() says so when the chain has not been post-processed", {
  out <- utils::capture.output(print(read_geneland(gl_fixture_run_raw())))
  expect_true(any(grepl("not run", out)))
})

test_that("a missing output file gives an informative error, not a low-level one", {
  ## Copy a run and delete one output to simulate write.* = FALSE.
  src <- gl_fixture_run_raw()
  dst <- gl_new_dir("missing-file")
  file.copy(list.files(src, full.names = TRUE), dst)
  file.remove(file.path(dst, "drifts.txt"))

  run <- read_geneland(dst)
  err <- tryCatch(autoplot(run, "drift"), error = conditionMessage)
  expect_match(err, "drifts.txt")
  expect_match(err, "write\\.\\* option", perl = FALSE)
})
