test_that("FormatGenotypes() recodes alleles consecutively from one", {
  dat <- gl_example_data()
  fg <- FormatGenotypes(genotypes = dat$geno, ploidy = 2)

  expect_named(fg, c("genotypes", "allele.numbers"))
  expect_equal(dim(fg$genotypes), dim(dat$geno))
  expect_length(fg$allele.numbers, ncol(dat$geno) / 2)

  ## Every locus must use 1..k with no gaps, which is what the sampler assumes.
  for (loc in seq_len(ncol(dat$geno) / 2)) {
    a <- as.vector(fg$genotypes[, (2 * loc - 1):(2 * loc)])
    a <- a[!is.na(a)]
    expect_equal(sort(unique(a)), seq_len(fg$allele.numbers[loc]),
                 info = paste("locus", loc))
  }
})

test_that("FormatGenotypes() preserves which entries are missing", {
  dat <- gl_example_data()
  geno <- dat$geno
  geno[1, 1] <- NA
  geno[5, 4] <- NA
  fg <- FormatGenotypes(genotypes = geno, ploidy = 2)
  expect_equal(is.na(fg$genotypes), is.na(geno))
})

test_that("gl2gp() writes a Genepop file with one line per individual", {
  dat <- gl_example_data()
  f <- tempfile(fileext = ".txt")
  gl2gp(coordinates = dat$coord, genotypes = dat$geno, file = f)

  expect_true(file.exists(f))
  lines <- readLines(f, warn = FALSE)
  nloc <- ncol(dat$geno) / 2
  nindiv <- nrow(dat$geno)

  ## gl2gp() opens a new Genepop "Pop" section whenever the coordinates change,
  ## so georeferenced individuals at distinct sites each get their own section:
  ## header + one line per locus + (Pop + individual) per individual.
  npop_lines <- sum(trimws(lines) == "Pop")
  expect_equal(npop_lines, nindiv)
  expect_equal(length(lines), 1 + nloc + 2 * nindiv)
  expect_equal(lines[1], "Header line")
})

test_that("gl2gp() codes missing genotypes as zeros", {
  dat <- gl_example_data()
  geno <- dat$geno
  geno[1, 1:2] <- NA
  f <- tempfile(fileext = ".txt")
  gl2gp(coordinates = dat$coord, genotypes = geno, file = f)

  lines <- readLines(f, warn = FALSE)
  first_indiv <- lines[1 + ncol(geno) / 2 + 1 + 1]
  expect_match(first_indiv, "0")
})

test_that("nullify() introduces null alleles without changing the shape", {
  dat <- gl_example_data()
  set.seed(3)
  res <- nullify(genotypes = dat$geno, nall.null = 1, nloc.null = 2)
  g <- if (is.list(res)) res$genotypes else res
  expect_equal(dim(g), dim(dat$geno))
})
