## Re-exported so that autoplot(run) works after library(Geneland) alone,
## without also requiring library(ggplot2).

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot
