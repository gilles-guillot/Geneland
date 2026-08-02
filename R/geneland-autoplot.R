## ggplot2 rendering of a Geneland run ------------------------------------

#' ggplot2 displays of a Geneland run
#'
#' Builds the standard Geneland diagnostics and maps as \pkg{ggplot2} objects.
#' Unlike the older `Plot*` functions these return a plot rather than drawing
#' one, so they can be modified, arranged and saved by the caller:
#'
#' ```r
#' p <- autoplot(run, "npop")
#' p + ggplot2::labs(title = "My dataset")
#' ggplot2::ggsave("npop.pdf", p, width = 6, height = 4)
#' ```
#'
#' They also leave the graphics device alone: no `dev.new()`, no `par()`
#' changes, and no second copy of the plotting code for writing to file.
#'
#' @param object A `geneland_run` object from [read_geneland()].
#' @param type Which display to produce:
#'   \describe{
#'     \item{`"npop"`}{trace of the number of clusters along the chain.}
#'     \item{`"npop_post"`}{posterior distribution of the number of clusters,
#'       after burn-in, with the modal value highlighted.}
#'     \item{`"ntile"`}{trace of the number of tiles in the tessellation.}
#'     \item{`"ntile_post"`}{posterior distribution of the number of tiles.}
#'     \item{`"drift"`}{drift factor traces, one panel per cluster
#'       (`freq.model = "Correlated"` only).}
#'     \item{`"loglik"`}{log-likelihood trace.}
#'     \item{`"logpost"`}{log-posterior-density trace.}
#'     \item{`"rate"`}{trace of the Poisson process rate.}
#'     \item{`"map"`}{map of modal cluster membership. Requires
#'       [PostProcessChain()].}
#'     \item{`"proba"`}{posterior probability of membership, one panel per
#'       cluster. Requires [PostProcessChain()].}
#'   }
#' @param burnin Number of *saved* iterations to discard. On traces the
#'   discarded region is shaded rather than hidden, so that a burn-in which is
#'   too short stays visible. Defaults to the value used by
#'   [PostProcessChain()] when that is known, otherwise 0.
#' @param coordinates Optional two-column matrix of individual coordinates,
#'   overlaid as points on `"map"` and `"proba"`.
#' @param ... Ignored.
#'
#' @return A `ggplot` object.
#'
#' @seealso [read_geneland()], [theme_geneland()], [geneland_pal()].
#'   The base-graphics functions ([Plotnpop()], [PosteriorMode()] and the
#'   other `Plot*` functions) are unchanged and remain available.
#'
#' @examples
#' \dontrun{
#' run <- read_geneland("/path/to/my/run/")
#' autoplot(run, "npop")
#' autoplot(run, "map", coordinates = coord)
#'
#' # Several panels side by side, if patchwork is installed:
#' # (autoplot(run, "npop") | autoplot(run, "npop_post"))
#' }
#' @importFrom ggplot2 autoplot
#' @export
autoplot.geneland_run <- function(object,
                                  type = c("npop", "npop_post",
                                           "ntile", "ntile_post",
                                           "drift", "loglik", "logpost",
                                           "rate", "map", "proba"),
                                  burnin = NULL,
                                  coordinates = NULL,
                                  ...) {
  type <- match.arg(type)
  if (is.null(burnin))
    burnin <- gl_param(object, "burnin", 0, "postprocess")

  switch(type,
    npop       = gl_plot_trace(object, "populations.numbers.txt", burnin,
                               "Number of clusters",
                               "Number of clusters along the chain",
                               integer_y = TRUE),
    npop_post  = gl_plot_posterior(object, "populations.numbers.txt", burnin,
                                   "Number of clusters",
                                   "Posterior distribution of the number of clusters",
                                   support = seq.int(gl_param(object, "npopmin", 1),
                                                     gl_param(object, "npopmax", 1))),
    ntile      = gl_plot_trace(object, "nuclei.numbers.txt", burnin,
                               "Number of tiles",
                               "Number of tiles in the tessellation",
                               integer_y = TRUE),
    ntile_post = gl_plot_posterior(object, "nuclei.numbers.txt", burnin,
                                   "Number of tiles",
                                   "Posterior distribution of the number of tiles"),
    loglik     = gl_plot_trace(object, "log.likelihood.txt", burnin,
                               "Log-likelihood", "Log-likelihood along the chain"),
    logpost    = gl_plot_trace(object, "log.posterior.density.txt", burnin,
                               "Log-posterior density",
                               "Log-posterior density along the chain"),
    rate       = gl_plot_trace(object, "Poisson.process.rate.txt", burnin,
                               "Poisson process rate",
                               "Rate of the Poisson process along the chain"),
    drift      = gl_plot_drift(object, burnin),
    map        = gl_plot_map(object, coordinates),
    proba      = gl_plot_proba(object, coordinates)
  )
}

#' @rdname autoplot.geneland_run
#' @param x A `geneland_run` object.
#' @param y Ignored.
#' @export
plot.geneland_run <- function(x, y, ...) {
  print(autoplot(x, ...))
  invisible(x)
}

## Building blocks -------------------------------------------------------

## Shade the burn-in rather than silently dropping it: a burn-in that is too
## short is then visible in the plot instead of hidden by it. The region is
## explained in the caption rather than by an inline label, which would be
## repeated in every panel of a faceted plot.
gl_burnin_layer <- function(burnin, thinning) {
  if (burnin <= 0) return(NULL)
  ggplot2::annotate("rect",
                    xmin = -Inf, xmax = burnin * thinning,
                    ymin = -Inf, ymax = Inf,
                    fill = gl_ink$grid, alpha = 0.55)
}

gl_burnin_caption <- function(burnin, thinning) {
  if (burnin <= 0) return(NULL)
  sprintf("Shaded: burn-in of %d saved iterations (%d MCMC iterations).",
          burnin, burnin * thinning)
}

gl_plot_trace <- function(run, file, burnin, ylab, title, integer_y = FALSE) {
  df <- gl_chain_df(run, file)
  burnin <- gl_check_burnin(burnin, nrow(df))
  thin <- gl_param(run, "thinning", 1)

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$iteration, .data$value)) +
    gl_burnin_layer(burnin, thin) +
    ggplot2::geom_line(linewidth = 0.4, colour = gl_categorical[1]) +
    ggplot2::labs(x = "MCMC iteration", y = ylab, title = title,
                  caption = gl_burnin_caption(burnin, thin)) +
    theme_geneland()

  if (integer_y)
    p <- p + ggplot2::scale_y_continuous(breaks = gl_int_breaks)
  p
}

## Integer-only axis breaks: a "3.5th cluster" is not a thing.
gl_int_breaks <- function(lims) {
  b <- unique(floor(pretty(lims)))
  b[b == as.integer(b)]
}

gl_plot_posterior <- function(run, file, burnin, xlab, title, support = NULL) {
  df <- gl_chain_df(run, file)
  burnin <- gl_check_burnin(burnin, nrow(df))
  kept <- df$value[seq.int(burnin + 1L, nrow(df))]

  tab <- as.data.frame(table(value = kept), stringsAsFactors = FALSE)
  tab$value <- as.numeric(tab$value)
  tab$prob <- tab$Freq / sum(tab$Freq)

  ## For the number of clusters, show the whole prior support so that a
  ## degenerate posterior reads as "all the mass on one value" rather than as
  ## a single bar filling the panel with nothing to compare it against.
  if (!is.null(support)) {
    missing <- setdiff(support, tab$value)
    if (length(missing))
      tab <- rbind(tab, data.frame(value = missing, Freq = 0L, prob = 0))
    tab <- tab[order(tab$value), ]
  }

  ## Emphasis, not eight hues: the story is which value is modal.
  tab$modal <- tab$prob == max(tab$prob)

  ## Never a number on every bar. Label the mode always; label the rest only
  ## when there are few enough bars for the labels not to collide.
  lab <- tab[tab$modal | nrow(tab) <= 8L, , drop = FALSE]
  lab <- lab[lab$prob > 0, , drop = FALSE]

  ggplot2::ggplot(tab, ggplot2::aes(.data$value, .data$prob)) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$modal), width = 0.7) +
    ggplot2::geom_text(data = lab,
                       mapping = ggplot2::aes(label = formatC(.data$prob,
                                                              format = "f",
                                                              digits = 2)),
                       vjust = -0.6, size = 3, colour = gl_ink$secondary) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = gl_categorical[1],
                                          `FALSE` = gl_ink$axis),
                               guide = "none") +
    ggplot2::scale_x_continuous(breaks = gl_int_breaks) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(x = xlab, y = "Posterior probability", title = title,
                  caption = sprintf("After a burn-in of %d saved iterations (%d MCMC iterations).",
                                    burnin, burnin * gl_param(run, "thinning", 1))) +
    theme_geneland() +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
}

gl_plot_drift <- function(run, burnin) {
  m <- gl_read_matrix(run, "drifts.txt")
  thin <- gl_param(run, "thinning", 1)
  burnin <- gl_check_burnin(burnin, nrow(m))

  ## Clusters that are not in the model at a given iteration are written as
  ## the sentinel -999. Left as a number the trace plunges off the panel;
  ## as NA the line simply stops, which is what actually happened.
  m[m == -999] <- NA_real_

  ## One panel per cluster rather than one colour per cluster: with npopmax
  ## clusters the traces overlap badly and colour alone cannot separate them.
  df <- data.frame(
    iteration = rep(seq_len(nrow(m)) * thin, times = ncol(m)),
    value     = as.vector(m),
    cluster   = factor(rep(seq_len(ncol(m)), each = nrow(m)),
                       labels = paste("Cluster", seq_len(ncol(m))))
  )

  ggplot2::ggplot(df, ggplot2::aes(.data$iteration, .data$value)) +
    gl_burnin_layer(burnin, thin) +
    ggplot2::geom_line(linewidth = 0.4, colour = gl_categorical[1]) +
    ggplot2::facet_wrap(~ .data$cluster) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "MCMC iteration", y = "Drift factor",
                  title = "Drift factors along the chain",
                  caption = gl_burnin_caption(burnin, thin)) +
    theme_geneland()
}

## Read the post-processed grid, erroring helpfully if it is not there.
gl_need_postprocess <- function(run, what) {
  if (is.null(run$postprocess))
    stop("The ", what, " needs the output of PostProcessChain(), which has ",
         "not been run for '", run$path, "'.", call. = FALSE)
}

gl_coord_layer <- function(coordinates) {
  if (is.null(coordinates)) return(NULL)
  coordinates <- as.matrix(coordinates)
  if (ncol(coordinates) != 2L)
    stop("'coordinates' must have exactly two columns.", call. = FALSE)
  df <- data.frame(x = coordinates[, 1], y = coordinates[, 2])
  ## A thin surface-coloured ring keeps the points legible on any fill.
  list(
    ggplot2::geom_point(data = df, ggplot2::aes(.data$x, .data$y),
                        inherit.aes = FALSE, shape = 21, size = 1.1,
                        stroke = 0.4, fill = gl_ink$primary,
                        colour = gl_ink$surface)
  )
}

gl_plot_map <- function(run, coordinates) {
  gl_need_postprocess(run, "membership map")
  m <- gl_read_matrix(run, "modal.pop.txt")
  df <- data.frame(x = m[, 1], y = m[, 2],
                   cluster = factor(as.integer(m[, 3])))
  k <- nlevels(df$cluster)

  ## Colour cannot carry cluster identity on a map: any two clusters may end
  ## up adjacent, and beyond three slots the palette is no longer separable
  ## under colour-vision deficiency. So label each region as well.
  ## Anchor each label to the cell of its own cluster nearest that cluster's
  ## centroid. Independent medians of x and y can fall outside a non-convex
  ## region, putting the label on top of a different cluster.
  labs <- do.call(rbind, lapply(split(df, df$cluster), function(d) {
    cx <- mean(d$x); cy <- mean(d$y)
    i <- which.min((d$x - cx)^2 + (d$y - cy)^2)
    data.frame(x = d$x[i], y = d$y[i], cluster = d$cluster[1])
  }))

  ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, fill = .data$cluster)) +
    ggplot2::geom_raster() +
    gl_coord_layer(coordinates) +
    ggplot2::geom_label(data = labs,
                        mapping = ggplot2::aes(x = .data$x, y = .data$y,
                                               label = .data$cluster),
                        inherit.aes = FALSE,
                        size = 3.2, fontface = "bold",
                        colour = gl_ink$primary, fill = gl_ink$surface,
                        alpha = 0.85, linewidth = 0,
                        label.padding = ggplot2::unit(0.12, "lines")) +
    scale_fill_geneland_d(name = "Cluster") +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(x = "x coordinate", y = "y coordinate",
                  title = "Modal cluster membership",
                  subtitle = sprintf("%d clusters; points are sampled individuals", k)) +
    theme_geneland() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank())
}

gl_plot_proba <- function(run, coordinates) {
  gl_need_postprocess(run, "membership probability maps")
  m <- gl_read_matrix(run, "proba.pop.membership.txt")
  npop <- ncol(m) - 2L
  if (npop < 1L)
    stop("proba.pop.membership.txt has no cluster columns.", call. = FALSE)

  df <- data.frame(
    x       = rep(m[, 1], times = npop),
    y       = rep(m[, 2], times = npop),
    prob    = as.vector(m[, -(1:2), drop = FALSE]),
    cluster = factor(rep(seq_len(npop), each = nrow(m)),
                     labels = paste("Cluster", seq_len(npop)))
  )

  ggplot2::ggplot(df, ggplot2::aes(.data$x, .data$y, fill = .data$prob)) +
    ggplot2::geom_raster() +
    gl_coord_layer(coordinates) +
    ggplot2::facet_wrap(~ .data$cluster) +
    scale_fill_geneland_c(name = "Posterior\nprobability") +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(x = "x coordinate", y = "y coordinate",
                  title = "Posterior probability of cluster membership") +
    theme_geneland() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank())
}
