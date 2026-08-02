#' Launch the Geneland graphical interface
#'
#' Starts a Shiny application covering the standard Geneland workflow: load
#' data, choose the model, run the sampler, inspect the chain, and view and
#' export the maps.
#'
#' The application writes everything under [tempdir()], so it never touches
#' your filespace unless you use one of its download buttons. Every figure it
#' draws comes from [autoplot.geneland_run()], i.e. the same code the console
#' API uses.
#'
#' The sampler runs inside the R process serving the app, so the interface is
#' unresponsive while a chain is running. For long chains, call [MCMC()] from
#' the console and then point the app (or [read_geneland()]) at the output
#' directory.
#'
#' @param launch.browser Passed to [shiny::runApp()]. `TRUE` (the default in
#'   an interactive session) opens the app in your web browser.
#' @param port Port to serve on. `NULL` lets Shiny choose.
#' @param ... Further arguments passed to [shiny::runApp()].
#'
#' @return Invoked for its side effect; does not return until the app is
#'   closed.
#'
#' @seealso [Geneland.GUI.tcltk()] for the legacy Tk interface.
#'
#' @examples
#' \dontrun{
#' run_geneland_app()
#' }
#' @export
run_geneland_app <- function(launch.browser = interactive(),
                             port = NULL, ...) {
  missing <- c("shiny", "bslib", "ggplot2")[
    !vapply(c("shiny", "bslib", "ggplot2"),
            requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("The graphical interface needs the following package(s): ",
         paste(missing, collapse = ", "), ".\n",
         'Install them with: install.packages(c("',
         paste(missing, collapse = '", "'), '"))', call. = FALSE)

  appdir <- system.file("shiny", package = "Geneland")
  if (!nzchar(appdir) || !file.exists(file.path(appdir, "app.R")))
    stop("The Shiny application was not found in the installed package.",
         call. = FALSE)

  shiny::runApp(appdir, launch.browser = launch.browser, port = port, ...)
}

#' Menu-driven interface to Geneland
#'
#' Launches the Geneland graphical interface. Since version 5.0.0 this is the
#' Shiny application described in [run_geneland_app()]; the original Tk
#' interface is still available as [Geneland.GUI.tcltk()].
#'
#' @param lib.loc Ignored, retained for backward compatibility.
#' @param ... Passed to [run_geneland_app()].
#'
#' @return Invoked for its side effect.
#'
#' @examples
#' \dontrun{
#' Geneland.GUI()
#' }
#' @export
Geneland.GUI <- function(lib.loc = NULL, ...) {
  if (!is.null(lib.loc))
    message("'lib.loc' is ignored by the Shiny interface.")
  message("Starting the Geneland Shiny interface. ",
          "For the previous Tk interface, use Geneland.GUI.tcltk().")
  run_geneland_app(...)
}
