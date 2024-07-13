#' Launch Shiny App for likelihoodAsy rstar Analysis
#'
#' This function launches a Shiny application that facilitates the setup and execution
#' of likelihoodAsy rstar analysis. The app allows users to upload a dataset, specify a
#' model and parameters of interest, and perform the analysis with the option to compute
#' confidence intervals for r* statistics.
#'
#' @return A Shiny app object that can be run locally.
#'
#' @examples
#' if (interactive()) {
#'   LA_rstar()
#' }
#'
#' @references
#' Pierce, D. A., & Bellio, R. (2017). Modern Likelihood-Frequentist Inference.
#' International Statistical Review / Revue Internationale de Statistique, 85(3),
#' 519–541. <doi:10.1111/insr.12232>
#'
#' Bellio R, Pierce D (2020). likelihoodAsy: Functions for Likelihood Asymptotics.
#' R package version 0.51, \url{https://CRAN.R-project.org/package=likelihoodAsy}.
#'
#' @export
LA_rstar <- function() {
  # UI
  ui <- shiny::fluidPage(
    theme = shinythemes::shinytheme("lumen"),
    shiny::titlePanel("Set Up likelihoodAsy rstar Analysis"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::fileInput("datafile", "Choose CSV File",
                         multiple = FALSE,
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
        shiny::uiOutput("rstar_input_ui"),
        shiny::actionButton("show_citations", "Citations")
      ),
      shiny::mainPanel(
        shiny::uiOutput("variables_title"),  # Placeholder for the title
        DT::dataTableOutput("variables_table"),
        shiny::uiOutput("fitmod_summary_header"),
        shiny::verbatimTextOutput("fitmod_summary_output"),
        shiny::uiOutput("rstar_summary_header"),
        shiny::verbatimTextOutput("rstar_summary_output"),
        # Conditional UI for rstarci summary
        shiny::conditionalPanel(
          condition = "input.rstar_ci == true",
          shiny::uiOutput("rstarci_summary_header"),
          shiny::verbatimTextOutput("rstarci_summary_output")
        ),
        shiny::uiOutput("citation_header"),
        shiny::verbatimTextOutput("citations_output")
      )
    )
  )

  # Server
  server <- function(input, output, session) {

    # reactive: Read the uploaded CSV file
    uploaded_data <- shiny::reactiveVal()
    shiny::observe({
      inFile <- input$datafile
      if (!is.null(inFile)) {
        data <- utils::read.csv(inFile$datapath, stringsAsFactors = TRUE)
        uploaded_data(data)
      }
    })

    # set available variables title
    output$variables_title <- shiny::renderUI({
      if (!is.null(uploaded_data()) && nrow(uploaded_data()) > 0) {
        shiny::tags$h2("Available Variables")
      }
    })

    # create variables table
    output$variables_table <- DT::renderDataTable({
      shiny::req(uploaded_data())
      data <- uploaded_data()
      df <- data.frame(Variable = names(data), Type = sapply(data, class))
      DT::datatable(df, editable = 'cell', options = list(pageLength = 5),
                    rownames = FALSE)
    })

    # handle variable type edits
    shiny::observeEvent(input$variables_table_cell_edit, {
      info <- input$variables_table_cell_edit
      shiny::req(uploaded_data())
      data <- uploaded_data()
      row_number <- info$row
      new_value <- info$value

      if (info$col == 0) {
        tryCatch({
          names(data)[row_number] <- new_value
          uploaded_data(data)
        }, error = function(e) {
          shiny::showNotification(
            paste("Error in changing variable name:", e$message),
            type = "error",
            duration = NULL
          )
        })
      }

      if (info$col == 1) {
        variable_name <- names(data)[row_number]
        tryCatch({
          if (new_value == "factor") {
            data[[variable_name]] <- as.factor(data[[variable_name]])
          } else if (new_value == "numeric") {
            data[[variable_name]] <- as.numeric(data[[variable_name]])
          } else if (new_value == "integer") {
            data[[variable_name]] <- as.integer(data[[variable_name]])
          } else if (new_value == "double") {
            data[[variable_name]] <- as.double(data[[variable_name]])
          } else if (new_value == "character") {
            data[[variable_name]] <- as.character(data[[variable_name]])
          } else {
            stop("New data type must be one of the following: factor, numeric, integer, double, character")
          }
          uploaded_data(data)
        }, error = function(e) {
          shiny::showNotification(
            paste("Error in changing data type:", e$message),
            type = "error",
            duration = NULL
          )
        })
      }
    })

    # setup rstar ui
    shiny::observe({
      if (!is.null(uploaded_data())) {
        output$rstar_input_ui <- shiny::renderUI({
          shiny::tagList(
            shiny::tags$div(title = "Select the type of model",
                            shiny::selectInput("model", "Model",
                                               choices = c("logistic", "linear", "poisson"), selected = "linear")),
            shiny::tags$div(title = "Enter a description of the model to be fitted",
                            shiny::textInput("formula", "Formula", value = "")),
            shiny::tags$div(title = "Describe the parameter of interest",
                            shiny::textInput("psidesc", "Psi Description", value = "Coefficient of Interest")),
            shiny::tags$div(title = "Specify the value of the parameter of interest under testing",
                            shiny::numericInput("psival", "Psi Value", value = 0)),
            shiny::tags$div(title = "Specify the index of the parameter of interest in the theta vector",
                            shiny::numericInput("fpsi", "Fpsi", value = 2)),
            shiny::tags$div(title = "Should the confidence interval for r* be calculated?",
                            shiny::checkboxInput("rstar_ci", "Calculate r* CI", value = FALSE)),
            # Additional parameters as needed
            shiny::tags$div(title = "Monte Carlo replicates used for computing the r* statistic",
                            shiny::numericInput("R", "Monte Carlo Replicates (R)", value = 1000)),
            shiny::tags$div(title = "Optional random seed for the Monte Carlo computation",
                            shiny::numericInput("seed", "Random Seed", value = NA)),
            shiny::tags$div(title = "Should computation trace be printed?",
                            shiny::checkboxInput("trace", "Print Computation Trace", value = FALSE)),
            shiny::tags$div(title = "Skip computation of the r* statistic",
                            shiny::checkboxInput("ronly", "Skip r* Computation", value = FALSE)),
            shiny::tags$div(title = "Constrained optimizer for maximizing the log likelihood function",
                            shiny::selectInput("constr_opt", "Constrained Optimizer",
                                               choices = c("solnp", "alabama"), selected = "solnp")),
            shiny::actionButton("run_analysis", "Run Analysis")
          )
        })
      }
    })

    # set up reactive values for BF analysis
    rstar_analysis_done <- shiny::reactiveVal(FALSE)
    rstar_result <- shiny::reactiveVal()

    # run BF_for_everyone analysis
    shiny::observeEvent(input$run_analysis, {
      shiny::req(uploaded_data(), input$formula, input$psival, input$fpsi)

      if (is.na(input$seed)) { set_seed <- NULL } else {set_seed <- input$seed}

      tryCatch({
        result <- rstar_glm(
          .model = input$model,
          .formula = stats::as.formula(input$formula),
          .data = uploaded_data(),
          .psidesc = input$psidesc,
          .psival = input$psival,
          .fpsi = input$fpsi,
          .rstar.ci = input$rstar_ci,
          R = input$R,
          seed = set_seed,
          trace = input$trace,
          ronly = input$ronly,
          constr.opt = input$constr_opt
        )

        # use reactive values to flag that analysis is done
        rstar_result(result)
        rstar_analysis_done(TRUE)

        # set output values related to BF analysis
        output$rstar_summary_output <- shiny::renderPrint({ summary(result$rs) })
        output$rstarci_summary_output <- shiny::renderPrint({ summary(result$rs_ci) })
        output$fitmod_summary_output <- shiny::renderPrint({ summary(result$fit_glm) })
      }, error = function(e) {
        shiny::showNotification(
          paste("Error:", e$message),
          type = "error",
          duration = NULL
        )
        rstar_analysis_done(FALSE)
      })
    })

    # setup rstar summary title
    output$rstar_summary_header <- shiny::renderUI({
      if (rstar_analysis_done()) {
        shiny::tags$h2("rstar Summary")
      }
    })

    # setup rstar.ci  title
    output$rstarci_summary_header <- shiny::renderUI({
      if (rstar_analysis_done()) {
        shiny::tags$h2("rstar CI Summary")
      }
    })

    # setup fitted model  title
    output$fitmod_summary_header <- shiny::renderUI({
      if (rstar_analysis_done()) {
        shiny::tags$h2("Fitted Model Summary")
      }
    })


    # Initialize citations_text as an empty string
    citations_text <- shiny::reactiveVal("")

    shiny::observeEvent(input$show_citations, {
      # Get the formatted citations
      likelihoodAsy_citation <- format_citation(utils::citation("likelihoodAsy"))
      holi_citation <- format_citation(utils::citation("holi"))

      citations <- paste(
        "Statistical Methods:",
        "Pierce, D.A. and Bellio, R. (2017). Modern likelihood-frequentist inference. International Statistical Review, 85, 519-541.",
        "",
        "likelihoodAsy Package:",
        likelihoodAsy_citation,
        "",
        "Web Application:",
        holi_citation,
        sep = "\n"
      )
      citations_text(citations)
    })


    # Render the citations output
    output$citations_output <- shiny::renderText({
      citations_text()
    })

    output$citation_header <- shiny::renderUI({
      shiny::req(citations_text())
        shiny::tags$h2("Citations")
    })

  }

  shiny::shinyApp(ui = ui, server = server)
}

#' Format Citation
#'
#' This internal function formats a citation object into a readable string.
#' The function extracts relevant information such as the title, author,
#' year, address, and URL from the citation object and formats it into a
#' standardized citation format.
#'
#' @param cit A citation object typically obtained from `citation()`.
#'
#' @return A character string with the formatted citation.
#'
#' @keywords internal
format_citation <- function(cit) {
  title <- cit$title
  author <- if (is.null(cit$author)) {
    cit$organization
  } else {
    paste(sapply(cit$author, function(a) paste(a$given, a$family)), collapse = ", ")
  }
  year <- cit$year
  address <- cit$address
  url <- cit$url

  formatted_cit <- paste0(
    author, " (", year, "). ",
    title, ". ",
    "Retrieved from ", url, ". ",
    address
  )

  formatted_cit
}
