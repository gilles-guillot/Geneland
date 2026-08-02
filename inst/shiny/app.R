## Geneland Shiny interface
##
## Launched by Geneland::run_geneland_app(). Kept as a single self-contained
## file so it can be read end to end.
##
## Design notes:
##  * Every plot is produced by autoplot() on a geneland_run object, i.e. the
##    same code path as the console API. There is deliberately no second
##    plotting implementation here.
##  * All output is written under tempdir(); nothing touches the user's
##    filespace unless they explicitly download it.

library(shiny)
library(bslib)
library(Geneland)

MAX_UPLOAD_MB <- 64
options(shiny.maxRequestSize = MAX_UPLOAD_MB * 1024^2)

## ---------------------------------------------------------------- helpers

## Read a whitespace/tab/comma-delimited numeric matrix, tolerating either.
read_matrix_file <- function(path, header) {
  first <- readLines(path, n = 1L, warn = FALSE)
  sep <- if (grepl(",", first)) "," else ""
  as.matrix(utils::read.table(path, header = header, sep = sep))
}

example_path <- function(f) system.file("extdata", f, package = "Geneland")

## ------------------------------------------------------------------- ui

ui <- page_navbar(
  title = "Geneland",
  theme = bs_theme(version = 5, preset = "flatly"),
  id = "nav",

  ## ---- 1. Data -------------------------------------------------------
  nav_panel(
    "1. Data",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        h5("Input data"),
        radioButtons("src", NULL,
                     c("Example dataset" = "example",
                       "Upload my own"   = "upload")),
        conditionalPanel(
          "input.src == 'upload'",
          fileInput("f_coord", "Coordinates (2 columns)",
                    accept = c(".txt", ".csv", ".dat")),
          fileInput("f_geno", "Genotypes (2 columns per locus)",
                    accept = c(".txt", ".csv", ".dat")),
          checkboxInput("header", "Files have a header row", FALSE)
        ),
        actionButton("load", "Load data", class = "btn-primary w-100"),
        hr(),
        helpText(sprintf("Uploads are capped at %d MB and are written to a ",
                         MAX_UPLOAD_MB),
                 "temporary directory that is removed when the session ends.")
      ),
      card(card_header("Summary"), verbatimTextOutput("data_summary")),
      layout_columns(
        card(card_header("Sampling locations"), plotOutput("data_map", height = 340)),
        card(card_header("Genotypes (first rows)"), tableOutput("geno_head"))
      )
    )
  ),

  ## ---- 2. Model ------------------------------------------------------
  nav_panel(
    "2. Model",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        h5("Model"),
        selectInput("freq_model", "Allele frequency model",
                    c("Uncorrelated", "Correlated")),
        checkboxInput("spatial", "Spatial (Poisson-Voronoi prior)", TRUE),
        checkboxInput("varnpop", "Treat the number of clusters as unknown", TRUE),
        checkboxInput("filterna", "Filter null alleles", FALSE),
        hr(),
        h5("Number of clusters"),
        numericInput("npopmin", "Minimum", 1, min = 1, step = 1),
        numericInput("npopinit", "Initial", 5, min = 1, step = 1),
        numericInput("npopmax", "Maximum", 8, min = 1, step = 1),
        hr(),
        h5("Chain"),
        numericInput("nit", "Iterations", 20000, min = 100, step = 1000),
        numericInput("thinning", "Thinning", 100, min = 1, step = 10),
        numericInput("delta_coord", "Coordinate uncertainty", 0, min = 0)
      ),
      card(
        card_header("What these settings mean"),
        htmlOutput("model_help")
      ),
      card(card_header("Derived quantities"), verbatimTextOutput("derived"))
    )
  ),

  ## ---- 3. Run --------------------------------------------------------
  nav_panel(
    "3. Run",
    layout_sidebar(
      sidebar = sidebar(
        width = 340,
        h5("Run the sampler"),
        actionButton("run", "Run MCMC", class = "btn-success w-100"),
        hr(),
        h5("Post-processing"),
        numericInput("burnin", "Burn-in (saved iterations)", 100, min = 0),
        numericInput("nxdom", "Grid columns", 100, min = 10, step = 10),
        numericInput("nydom", "Grid rows", 100, min = 10, step = 10),
        actionButton("post", "Post-process", class = "btn-primary w-100"),
        hr(),
        helpText("The sampler runs in this R process, so the interface is ",
                 "unresponsive while it works. Long chains are better ",
                 "launched from the console with MCMC().")
      ),
      card(card_header("Status"), verbatimTextOutput("run_status")),
      card(card_header("Run"), verbatimTextOutput("run_info"))
    )
  ),

  ## ---- 4. Diagnostics ------------------------------------------------
  nav_panel(
    "4. Diagnostics",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        selectInput("diag", "Display",
                    c("Number of clusters (trace)"     = "npop",
                      "Number of clusters (posterior)" = "npop_post",
                      "Number of tiles (trace)"        = "ntile",
                      "Number of tiles (posterior)"    = "ntile_post",
                      "Drift factors"                  = "drift",
                      "Log-likelihood"                 = "loglik",
                      "Log-posterior density"          = "logpost",
                      "Poisson process rate"           = "rate")),
        numericInput("diag_burnin", "Burn-in (saved iterations)", 100, min = 0),
        downloadButton("dl_diag", "Download PDF", class = "w-100")
      ),
      card(full_screen = TRUE, plotOutput("diag_plot", height = 560))
    )
  ),

  ## ---- 5. Maps -------------------------------------------------------
  nav_panel(
    "5. Maps",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        radioButtons("maptype", "Display",
                     c("Modal cluster membership" = "map",
                       "Membership probabilities" = "proba")),
        checkboxInput("show_pts", "Show sampled individuals", TRUE),
        downloadButton("dl_map", "Download PDF", class = "w-100"),
        hr(),
        helpText("Cluster numbers are drawn on the map as well as colours, ",
                 "so the clusters remain distinguishable with colour-vision ",
                 "deficiency.")
      ),
      card(full_screen = TRUE, plotOutput("map_plot", height = 620))
    )
  ),

  ## ---- 6. Results ----------------------------------------------------
  nav_panel(
    "6. Results",
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        h5("Export"),
        downloadButton("dl_all", "Download run directory (.zip)",
                       class = "btn-primary w-100"),
        hr(),
        helpText("The archive holds every file written by MCMC() and ",
                 "PostProcessChain(), so the run can be reopened later with ",
                 "read_geneland().")
      ),
      card(card_header("Estimated cluster membership"),
           tableOutput("membership")),
      card(card_header("Files written"), verbatimTextOutput("files"))
    )
  ),

  nav_spacer(),
  nav_item(tags$span(class = "navbar-text small",
                     textOutput("version", inline = TRUE)))
)

## --------------------------------------------------------------- server

server <- function(input, output, session) {

  rv <- reactiveValues(
    coord = NULL, geno = NULL, path = NULL, run = NULL,
    status = "No data loaded.", posted = FALSE
  )

  output$version <- renderText(
    paste("Geneland", as.character(utils::packageVersion("Geneland"))))

  ## ---- data ----------------------------------------------------------
  observeEvent(input$load, {
    tryCatch({
      ## Validate into locals and only publish to rv once everything passes:
      ## a rejected file must not become the active dataset.
      if (identical(input$src, "example")) {
        coord <- as.matrix(utils::read.table(example_path("coordinates.txt")))
        geno  <- as.matrix(utils::read.table(example_path("genotypes.txt")))
      } else {
        req(input$f_coord, input$f_geno)
        coord <- read_matrix_file(input$f_coord$datapath, input$header)
        geno  <- read_matrix_file(input$f_geno$datapath, input$header)
      }
      if (ncol(coord) != 2L)
        stop("The coordinates file must have exactly 2 columns; it has ",
             ncol(coord), ".")
      if (nrow(coord) != nrow(geno))
        stop("Coordinates have ", nrow(coord), " rows but genotypes have ",
             nrow(geno), ". They must match, one row per individual.")
      if (ncol(geno) %% 2L != 0L)
        stop("Diploid codominant genotypes need 2 columns per locus; the file ",
             "has an odd number of columns (", ncol(geno), ").")

      rv$coord <- coord; rv$geno <- geno
      rv$path <- NULL; rv$run <- NULL; rv$posted <- FALSE
      rv$status <- "Data loaded. Set the model, then run the sampler."
      updateNumericInput(session, "npopinit", value = 5)
      showNotification(sprintf("Loaded %d individuals, %d loci.",
                               nrow(rv$geno), ncol(rv$geno) / 2),
                       type = "message")
    }, error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL)
    })
  })

  output$data_summary <- renderPrint({
    if (is.null(rv$coord)) return(cat("No data loaded yet.\n"))
    cat("Individuals :", nrow(rv$coord), "\n")
    cat("Loci        :", ncol(rv$geno) / 2, "\n")
    cat("Missing     :", sum(is.na(rv$geno)), "genotype entries\n")
    cat("x range     :", paste(signif(range(rv$coord[, 1]), 4), collapse = " to "), "\n")
    cat("y range     :", paste(signif(range(rv$coord[, 2]), 4), collapse = " to "), "\n")
  })

  output$geno_head <- renderTable({
    if (is.null(rv$geno)) return(NULL)
    utils::head(as.data.frame(rv$geno[, seq_len(min(10, ncol(rv$geno))), drop = FALSE]), 6)
  }, rownames = TRUE)

  ## A blank panel reads as a broken app; say what the user should do instead.
  placeholder <- function(msg) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = msg, colour = "#898781") +
      ggplot2::theme_void()
  }

  output$data_map <- renderPlot({
    if (is.null(rv$coord)) return(placeholder("Load data to see the sampling locations."))
    df <- data.frame(x = rv$coord[, 1], y = rv$coord[, 2])
    ggplot2::ggplot(df, ggplot2::aes(x = .data$x, y = .data$y)) +
      ggplot2::geom_point(size = 1.8, colour = "#2a78d6") +
      ggplot2::coord_equal() +
      ggplot2::labs(x = "x coordinate", y = "y coordinate") +
      Geneland::theme_geneland()
  })

  ## ---- model ---------------------------------------------------------
  output$model_help <- renderUI({
    HTML(paste0(
      "<p><b>Uncorrelated</b> gives each cluster an independent prior on ",
      "allele frequencies. It is the safer default and the only model that ",
      "supports null-allele filtering.</p>",
      "<p><b>Correlated</b> draws frequencies around a common ancestral ",
      "population. It has more power to detect weak structure, but is more ",
      "prone to reporting spurious clusters when its assumptions fail.</p>",
      "<p>Choose <b>maximum</b> larger than any number of clusters you could ",
      "reasonably expect. If the trace on the Diagnostics tab presses against ",
      "that ceiling, raise it and rerun.</p>"))
  })

  derived <- reactive({
    req(rv$coord)
    n <- nrow(rv$coord)
    list(nindiv = n, rate_max = n, nb_nuclei_max = 3 * n,
         nsaved = floor(input$nit / max(1, input$thinning)))
  })

  output$derived <- renderPrint({
    if (is.null(rv$coord)) return(cat("Load data on tab 1 first.\n"))
    d <- derived()
    cat("Poisson rate max   :", d$rate_max, " (= number of individuals)\n")
    cat("Max nuclei         :", d$nb_nuclei_max, " (= 3 x rate max)\n")
    cat("Saved iterations   :", d$nsaved, "\n")
    if (d$nsaved < 100)
      cat("\nNote: fewer than 100 saved states makes the posterior summaries noisy.\n")
  })

  ## ---- run -----------------------------------------------------------
  observeEvent(input$run, {
    if (is.null(rv$coord)) {
      showNotification("Load data first.", type = "warning"); return()
    }
    if (input$npopmin > input$npopinit || input$npopinit > input$npopmax) {
      showNotification("Require minimum <= initial <= maximum.", type = "error")
      return()
    }
    d <- derived()
    path <- file.path(tempdir(), format(Sys.time(), "geneland-%Y%m%d-%H%M%S"))
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path <- paste0(path, .Platform$file.sep)

    withProgress(message = "Running MCMC", value = 0.1, {
      res <- try(utils::capture.output(
        Geneland::MCMC(
          coordinates = rv$coord, geno.dip.codom = rv$geno, path.mcmc = path,
          rate.max = d$rate_max, delta.coord = input$delta_coord,
          npopmin = input$npopmin, npopinit = input$npopinit,
          npopmax = input$npopmax, nb.nuclei.max = d$nb_nuclei_max,
          nit = input$nit, thinning = input$thinning,
          freq.model = input$freq_model, varnpop = input$varnpop,
          spatial = input$spatial, filter.null.alleles = input$filterna)),
        silent = TRUE)
      setProgress(1)
    })

    if (inherits(res, "try-error")) {
      rv$status <- paste("MCMC failed:", conditionMessage(attr(res, "condition")))
      showNotification(rv$status, type = "error", duration = NULL)
      return()
    }
    rv$path <- path
    rv$run <- Geneland::read_geneland(path)
    rv$posted <- FALSE
    rv$status <- "MCMC finished. Post-process to obtain the maps."
    updateNumericInput(session, "burnin", value = max(1, floor(d$nsaved / 5)))
    updateNumericInput(session, "diag_burnin", value = max(1, floor(d$nsaved / 5)))
    showNotification("MCMC finished.", type = "message")
    nav_select("nav", "4. Diagnostics")
  })

  observeEvent(input$post, {
    if (is.null(rv$path)) {
      showNotification("Run the sampler first.", type = "warning"); return()
    }
    withProgress(message = "Post-processing", value = 0.1, {
      res <- try(utils::capture.output(
        Geneland::PostProcessChain(coordinates = rv$coord, path.mcmc = rv$path,
                                   nxdom = input$nxdom, nydom = input$nydom,
                                   burnin = input$burnin)),
        silent = TRUE)
      setProgress(1)
    })
    if (inherits(res, "try-error")) {
      rv$status <- paste("Post-processing failed:",
                         conditionMessage(attr(res, "condition")))
      showNotification(rv$status, type = "error", duration = NULL)
      return()
    }
    rv$run <- Geneland::read_geneland(rv$path)
    rv$posted <- TRUE
    rv$status <- "Post-processing finished. The maps are ready."
    showNotification("Post-processing finished.", type = "message")
    nav_select("nav", "5. Maps")
  })

  output$run_status <- renderPrint(cat(rv$status, "\n"))

  output$run_info <- renderPrint({
    if (is.null(rv$run)) return(cat("Nothing has been run in this session.\n"))
    print(rv$run)
  })

  ## ---- plots ---------------------------------------------------------
  diag_plot <- reactive({
    req(rv$run)
    Geneland::autoplot(rv$run, input$diag, burnin = input$diag_burnin)
  })

  output$diag_plot <- renderPlot({
    if (is.null(rv$run)) return(placeholder("Run the sampler first (tab 3)."))
    diag_plot()
  }, res = 100)

  map_plot <- reactive({
    req(rv$run, rv$posted)
    Geneland::autoplot(rv$run, input$maptype,
                       coordinates = if (input$show_pts) rv$coord else NULL)
  })

  output$map_plot <- renderPlot({
    if (is.null(rv$run) || !rv$posted)
      return(placeholder("Post-process the chain first (tab 3)."))
    map_plot()
  }, res = 100)

  output$dl_diag <- downloadHandler(
    filename = function() paste0("geneland-", input$diag, ".pdf"),
    content = function(file)
      ggplot2::ggsave(file, diag_plot(), width = 8, height = 5.5))

  output$dl_map <- downloadHandler(
    filename = function() paste0("geneland-", input$maptype, ".pdf"),
    content = function(file)
      ggplot2::ggsave(file, map_plot(), width = 8, height = 6.5))

  ## ---- results -------------------------------------------------------
  output$membership <- renderTable({
    req(rv$run, rv$posted)
    f <- file.path(rv$run$path, "proba.pop.membership.indiv.txt")
    if (!file.exists(f)) return(NULL)
    m <- as.matrix(utils::read.table(f))
    p <- m[, -(1:2), drop = FALSE]
    data.frame(individual = seq_len(nrow(p)),
               x = m[, 1], y = m[, 2],
               cluster = max.col(p),
               probability = round(apply(p, 1, max), 3))
  }, digits = 3)

  output$files <- renderPrint({
    if (is.null(rv$path)) return(cat("Nothing written yet.\n"))
    cat(rv$path, "\n\n"); print(sort(basename(list.files(rv$path))))
  })

  output$dl_all <- downloadHandler(
    filename = function() "geneland-run.zip",
    content = function(file) {
      req(rv$path)
      utils::zip(file, list.files(rv$path, full.names = TRUE), flags = "-j9X")
    },
    contentType = "application/zip")
}

shinyApp(ui, server)
