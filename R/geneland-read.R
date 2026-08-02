## Reading a Geneland run directory into R objects -----------------------
##
## The historical plotting functions each re-open and re-parse the MCMC output
## files themselves. read_geneland() does that once and returns an object, so
## the ggplot2 layer in geneland-autoplot.R has a single, tested entry point.

#' Read a Geneland MCMC run
#'
#' Reads the output files written by [MCMC()] (and, if it has been run,
#' [PostProcessChain()]) into a single object that the plotting methods can
#' consume.
#'
#' The historical `Plot*` functions each re-read these files on every call.
#' Reading them once means the run metadata is validated in one place and the
#' plotting code does not need to know the file layout.
#'
#' @param path Path to the directory containing the output of [MCMC()], i.e.
#'   the same value passed as `path.mcmc`.
#'
#' @return An object of class `geneland_run`: a list with elements
#'   \describe{
#'     \item{`path`}{normalised path to the run directory.}
#'     \item{`params`}{named list of the run parameters, read from
#'       \file{parameters.txt} and coerced to numeric/logical where possible.}
#'     \item{`postprocess`}{named list of post-processing parameters read from
#'       \file{postprocess.parameters.txt}, or `NULL` if
#'       [PostProcessChain()] has not been run.}
#'   }
#'
#' @seealso [autoplot.geneland_run()] for the plots.
#' @examples
#' \dontrun{
#' run <- read_geneland("/path/to/my/run/")
#' run
#' }
#' @export
read_geneland <- function(path) {
  if (!is.character(path) || length(path) != 1L)
    stop("'path' must be a single character string.", call. = FALSE)
  if (!dir.exists(path))
    stop("Directory not found: ", path, call. = FALSE)

  pfile <- file.path(path, "parameters.txt")
  if (!file.exists(pfile))
    stop("'", path, "' does not look like a Geneland run directory: ",
         "parameters.txt is missing.", call. = FALSE)

  params <- gl_read_keyvalue(pfile)

  ppfile <- file.path(path, "postprocess.parameters.txt")
  postprocess <- if (file.exists(ppfile)) gl_read_keyvalue(ppfile) else NULL

  structure(
    list(path = normalizePath(path, mustWork = TRUE),
         params = params,
         postprocess = postprocess),
    class = "geneland_run")
}

#' @param x A `geneland_run` object.
#' @param ... Ignored.
#' @rdname read_geneland
#' @export
print.geneland_run <- function(x, ...) {
  cat("<geneland_run>\n")
  cat("  path       : ", x$path, "\n", sep = "")
  cat("  individuals: ", gl_param(x, "nindiv", NA), "\n", sep = "")
  cat("  iterations : ", gl_param(x, "nit", NA),
      "  (thinning ", gl_param(x, "thinning", NA),
      " -> ", gl_nsaved(x), " saved)\n", sep = "")
  cat("  npopmax    : ", gl_param(x, "npopmax", NA), "\n", sep = "")
  cat("  freq.model : ", gl_param(x, "freq.model", NA), "\n", sep = "")
  if (is.null(x$postprocess)) {
    cat("  postprocess: not run (PostProcessChain() gives the maps)\n")
  } else {
    cat("  postprocess: ", gl_param(x, "nxdom", NA, "postprocess"), " x ",
        gl_param(x, "nydom", NA, "postprocess"), " grid\n", sep = "")
  }
  invisible(x)
}

## Internal helpers ------------------------------------------------------

## Parse the "name : value" files Geneland writes.
gl_read_keyvalue <- function(file) {
  tab <- utils::read.table(file, stringsAsFactors = FALSE, fill = TRUE)
  vals <- as.character(tab[[ncol(tab)]])
  out <- lapply(vals, function(v) {
    if (v %in% c("TRUE", "FALSE")) return(as.logical(v))
    num <- suppressWarnings(as.numeric(v))
    if (!is.na(num)) num else v
  })
  stats::setNames(out, as.character(tab[[1L]]))
}

gl_param <- function(run, name, default = NULL, where = "params") {
  v <- run[[where]][[name]]
  if (is.null(v)) default else v
}

## Number of saved (thinned) iterations.
gl_nsaved <- function(run) {
  nit <- gl_param(run, "nit", NA)
  thin <- gl_param(run, "thinning", 1)
  if (is.na(nit)) NA_integer_ else as.integer(nit / thin)
}

gl_file <- function(run, name) {
  f <- file.path(run$path, name)
  if (!file.exists(f))
    stop("File '", name, "' not found in ", run$path,
         ".\nIt is produced only when the corresponding write.* option of ",
         "MCMC() (or PostProcessChain()) has been used.", call. = FALSE)
  f
}

## A chain stored one value per line.
gl_read_chain <- function(run, name) {
  scan(gl_file(run, name), quiet = TRUE)
}

gl_read_matrix <- function(run, name) {
  as.matrix(utils::read.table(gl_file(run, name)))
}

## Validate and normalise a burnin given in *saved* iterations.
gl_check_burnin <- function(burnin, n) {
  if (!is.numeric(burnin) || length(burnin) != 1L || is.na(burnin) || burnin < 0)
    stop("'burnin' must be a single non-negative number.", call. = FALSE)
  if (burnin >= n)
    stop("'burnin' (", burnin, ") must be smaller than the number of saved ",
         "iterations (", n, ").", call. = FALSE)
  as.integer(burnin)
}

## Long data frame for a scalar chain, with the iteration index expressed in
## MCMC iterations (not saved states) so the x axis means what it says.
gl_chain_df <- function(run, name) {
  v <- gl_read_chain(run, name)
  thin <- gl_param(run, "thinning", 1)
  data.frame(saved = seq_along(v),
             iteration = seq_along(v) * thin,
             value = v)
}
