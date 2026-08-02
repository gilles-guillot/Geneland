## Palette and theme for the ggplot2 layer --------------------------------
##
## Colour choices here are not a matter of taste; they were validated for
## colour-vision deficiency before being adopted. See the notes on each object.

## Chart chrome.
gl_ink <- list(
  primary   = "#0b0b0b",
  secondary = "#52514e",
  muted     = "#898781",
  grid      = "#e1e0d9",
  axis      = "#c3c2b7",
  surface   = "#fcfcfb"
)

## Categorical hues, in this fixed order. The ORDER is the colour-vision
## safety mechanism, not decoration: it is the order that clears the adjacent
## -pair separation gates. Never cycle or generate a 9th hue.
gl_categorical <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100",
                    "#e87ba4", "#008300", "#4a3aa7", "#e34948")

## Single-hue sequential ramp (light -> dark) for magnitudes such as posterior
## probabilities. Deliberately one hue: a rainbow ramp invents structure that
## is not in the data.
gl_sequential <- c("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5",
                   "#256abf", "#184f95", "#0d366b")

#' Geneland colour palette
#'
#' The categorical palette used for cluster identity. The slot order is fixed
#' because it is what makes adjacent clusters distinguishable under the common
#' forms of colour blindness; requesting colours in a different order, or
#' beyond the eighth, defeats that.
#'
#' Note that on a map any two clusters can end up adjacent. Only the first
#' three slots are separable under every pairing, so [autoplot.geneland_run()]
#' always draws the cluster number on the map as well: identity never rests on
#' colour alone.
#'
#' @param n Number of colours required (1 to 8).
#' @return A character vector of `n` hex colours.
#' @examples
#' geneland_pal(4)
#' @export
geneland_pal <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1)
    stop("'n' must be a single positive number.", call. = FALSE)
  n <- as.integer(n)
  if (n > length(gl_categorical))
    stop("The categorical palette has only ", length(gl_categorical),
         " slots; ", n, " were requested.\n",
         "Beyond that, colours cannot be told apart reliably. Use the ",
         "'proba' plot (one panel per cluster) instead.", call. = FALSE)
  gl_categorical[seq_len(n)]
}

#' Geneland ggplot2 theme
#'
#' A light theme with a hairline grid and recessive axes, so the data are the
#' most prominent thing in the panel.
#'
#' @param base_size Base font size in points.
#' @return A ggplot2 theme object.
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_geneland()
#' }
#' @export
theme_geneland <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(colour = gl_ink$primary,
                                               face = "bold",
                                               size = ggplot2::rel(1.05)),
      plot.subtitle    = ggplot2::element_text(colour = gl_ink$secondary,
                                               size = ggplot2::rel(0.9)),
      plot.caption     = ggplot2::element_text(colour = gl_ink$muted,
                                               size = ggplot2::rel(0.8),
                                               hjust = 0),
      axis.title       = ggplot2::element_text(colour = gl_ink$secondary),
      axis.text        = ggplot2::element_text(colour = gl_ink$muted),
      panel.grid.major = ggplot2::element_line(colour = gl_ink$grid,
                                               linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = gl_ink$surface,
                                               colour = NA),
      plot.background  = ggplot2::element_rect(fill = gl_ink$surface,
                                               colour = NA),
      strip.text       = ggplot2::element_text(colour = gl_ink$primary,
                                               face = "bold",
                                               size = ggplot2::rel(0.9)),
      legend.title     = ggplot2::element_text(colour = gl_ink$secondary),
      legend.text      = ggplot2::element_text(colour = gl_ink$muted),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
}

#' Geneland colour and fill scales
#'
#' Discrete scales use the validated cluster palette; the continuous fill scale
#' is the single-hue sequential ramp used for posterior probabilities.
#'
#' @param ... Passed to the underlying ggplot2 scale.
#' @param limits Passed to the underlying ggplot2 scale.
#' @return A ggplot2 scale object.
#' @name geneland_scales
#' @examples
#' \dontrun{
#' library(ggplot2)
#' ggplot(iris, aes(Sepal.Length, Sepal.Width, colour = Species)) +
#'   geom_point() + scale_colour_geneland_d() + theme_geneland()
#' }
NULL

#' @rdname geneland_scales
#' @export
scale_colour_geneland_d <- function(...) {
  ggplot2::discrete_scale("colour", palette = function(n) geneland_pal(n), ...)
}

#' @rdname geneland_scales
#' @export
scale_fill_geneland_d <- function(...) {
  ggplot2::discrete_scale("fill", palette = function(n) geneland_pal(n), ...)
}

#' @rdname geneland_scales
#' @export
scale_fill_geneland_c <- function(..., limits = c(0, 1)) {
  ggplot2::scale_fill_gradientn(colours = gl_sequential, limits = limits, ...)
}
