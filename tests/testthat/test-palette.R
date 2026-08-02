test_that("geneland_pal() returns the fixed slots in order", {
  expect_length(geneland_pal(1), 1)
  expect_length(geneland_pal(8), 8)
  ## The order is the colour-vision-safety mechanism: slot n must not change
  ## when more slots are requested.
  expect_identical(geneland_pal(3), geneland_pal(8)[1:3])
  expect_true(all(grepl("^#[0-9a-fA-F]{6}$", geneland_pal(8))))
  expect_equal(anyDuplicated(geneland_pal(8)), 0L)
})

test_that("geneland_pal() refuses to invent a ninth hue", {
  expect_error(geneland_pal(9), "only 8")
  expect_error(geneland_pal(0), "positive")
  expect_error(geneland_pal(-1), "positive")
  expect_error(geneland_pal("three"), "positive")
})

test_that("theme_geneland() is a usable theme", {
  th <- theme_geneland()
  expect_s3_class(th, "theme")
  expect_no_error(
    ggplot2::ggplot_build(
      ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) +
        ggplot2::geom_point() + th))
})

test_that("the scales attach and behave", {
  df <- data.frame(x = 1:6, y = 1:6, g = letters[1:6], v = seq(0, 1, length = 6))

  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = g)) +
    ggplot2::geom_point() + scale_colour_geneland_d()
  expect_no_error(ggplot2::ggplot_build(p))

  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = g)) +
    ggplot2::geom_col() + scale_fill_geneland_d()
  expect_no_error(ggplot2::ggplot_build(p))

  ## The continuous fill is pinned to [0, 1]: these are probabilities, and a
  ## free scale would make two panels incomparable.
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = v)) +
    ggplot2::geom_raster() + scale_fill_geneland_c()
  expect_equal(p$scales$get_scales("fill")$limits, c(0, 1))
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("a discrete scale past eight levels errors rather than cycling", {
  df <- data.frame(x = 1:9, y = 1:9, g = letters[1:9])
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, colour = g)) +
    ggplot2::geom_point() + scale_colour_geneland_d()
  expect_error(ggplot2::ggplot_build(p))
})
