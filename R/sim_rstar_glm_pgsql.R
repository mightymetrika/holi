sim_rstar_glm_pgsql <- function(dbname, datatable, host, port, user, password) {

  # Helper function to
  text_to_vector <- function(text_input) {
    as.numeric(unlist(strsplit(text_input, ",")))
  }

  # Helper function to save data to the database
  saveData <- function(data) {
    # Connect to the database
    pool <- pool::dbPool(
      drv = RPostgres::Postgres(),
      dbname = dbname,
      host = host,
      user = user,
      password = password,
      port = port
    )

    # Close pool on stop
    shiny::onStop(function() {
      pool::poolClose(pool)
    })

    # Convert NA to NaN in the data frame
    data[is.na(data)] <- NaN

    # Loop through rows of data and save to databse
    lapply(1:nrow(data), function(i){

      # get row i of the data
      row_data <- data[i, ]

      # Construct the update query by looping over the data fields
      query <- sprintf(
        "INSERT INTO %s (%s) VALUES ('%s')",
        datatable,
        paste(names(row_data), collapse = ", "),
        paste(row_data, collapse = "', '")
      )

      # Execute the query
      tryCatch({
        pool::dbExecute(pool, query)
      }, error = function(e) {
        print(paste("Error inserting row", i, ":", e))
      })
    })


  }

  # Helper function to load data from database
  loadData <- function() {

    # Connect to the database
    pool <- pool::dbPool(
      drv = RPostgres::Postgres(),
      dbname = dbname,
      host = host,
      user = user,
      password = password,
      port = port
    )

    # Close pool on stop
    shiny::onStop(function() {
      pool::poolClose(pool)
    })

    # Construct the fetching query
    sql <- sprintf("SELECT * FROM %s", datatable)
    query <- pool::sqlInterpolate(pool, sql)

    # Submit the fetch query and disconnect
    pool::dbGetQuery(pool, query)

  }

  # Define the UI
  ui <- shiny::fluidPage(
    shiny::titlePanel("rstar GLM Simulation"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        shiny::numericInput("n_sims", "Number of iterations:", 10),
        shiny::numericInput("alpha_level", "Alpha-level", 0.05),
        shiny::numericInput("n_main", "Main group size:", 20),
        shiny::numericInput("n_covariates", "Number of covariates:", 3),
        shiny::textInput("true_coef_main", "True coefficients (main):", ".5, -0.3, 0.2"),
        shiny::numericInput("n_control", "Control group size:", 20),
        shiny::textInput("true_coef_control", "True coefficients (control):", ".5, -0.3, 0.2"),
        shiny::numericInput("treatment_effect", "Treatment effect:", 0.5),
        shiny::selectInput("model", "Model",
                           choices = c("logistic", "linear", "poisson"), selected = "logistic"),
        shiny::textInput("skewness_main", "Covariate skew (main):", "0.08, 0.08, 0.08"),
        shiny::textInput("skewness_control", "Covariate skew (control):", "0.02, 0.02, 0.02"),
        shiny::textInput("Sigma_main", "Covariance matrix (main):", "stats::rWishart(1, 10, Sigma = stats::toeplitz((3:1)/(3*5)))[,,1]"),
        shiny::textInput("Sigma_control", "Covariance matrix (control):", "stats::rWishart(1, 10, Sigma = stats::toeplitz((3:1)/(3*5)))[,,1]"),
        shiny::actionButton("runSim", "Run Simulation"),
        shiny::actionButton("submit", "Submit"),
        shiny::br(),  # Add a line break
        shiny::br(),  # Add a line break
        shiny::downloadButton("downloadBtn", "Download Data")
      ),
      shiny::mainPanel(
        # # Conditionally display the Simulation Results header and table
        # shiny::uiOutput("simulation_results_header"),
        # DT::DTOutput("resultsTable"),
        # shiny::br(),  # Add a line break
        # shiny::br(),  # Add a line break
        # # Add a header for the responses table
        # shiny::div(
        #   shiny::h4("All Responses"),
        #   DT::DTOutput("responses")
        shiny::uiOutput("simulation_results_header"),
        DT::DTOutput("resultsTable"),
        shiny::br(),
        shiny::plotOutput("estimatePlot"),
        shiny::br(),
        shiny::plotOutput("pValuePlot"),
        shiny::br(),
        shiny::div(
          shiny::h4("All Responses"),
          DT::DTOutput("responses")
        )
      )
    )
  )

  # Define the server logic
  server <- function(input, output, session) {

    # # Render the UI for parameters based on the selected cell block
    # output$paramsUI <- shiny::renderUI({
    #   getUIParams_pgsql(input$cellBlock)
    # })

    # Reactive value to store the results
    results <- shiny::reactiveVal(data.frame())     #For display
    results_plot <- shiny::reactiveVal(NULL)     #For plot
    results_exp <- shiny::reactiveVal(data.frame()) #For export

    # Load data from the database on app start
    output$responses <- DT::renderDT({
      loadData()
    }, options = list(pageLength = 5))

    # Observe event for the run simulation button
    shiny::observeEvent(input$runSim, {

      # Handle null values for text to vector
      vec_null <- function(par_input = "") {
        if (is.na(par_input) || par_input == "") {
          return(NULL)
        } else {
          return(text_to_vector(par_input))
        }
      }

      # Handle null values for sigma matrix
      sig_null <- function(par_input = "") {
        if (is.na(par_input) || par_input == "") {
          return(NULL)
        } else {
          return(eval(parse(text = par_input)))
        }
      }

      # Run simulation
      # Call the simulation function with both user-provided and default parameters
      simResults <- run_sim_rstar_glm(n_sims = input$n_sims,
                                      alpha_level = input$alpha_level,
                                      n_main = input$n_main,
                                      n_covariates = input$n_covariates,
                                      true_coef_main = vec_null(input$true_coef_main),
                                      n_control = input$n_control,
                                      true_coef_control = vec_null(input$true_coef_control),
                                      treatment_effect = input$treatment_effect,
                                      model = input$model,
                                      skewness_main = vec_null(input$skewness_main),
                                      skewness_control = vec_null(input$skewness_control),
                                      Sigma_main = sig_null(input$Sigma_main),
                                      Sigma_control = sig_null(input$Sigma_control)
                                      )

      # Update the results reactive values
      results(simResults$summary)
      results_plot(simResults$results)
      results_exp(simResults$summary)
    })

    #Output the results table
    output$resultsTable <- DT::renderDT({
      results()
    }, options = list(pageLength = 5))

    # Add reactive plot generation
    plots <- shiny::reactive({
      shiny::req(results_plot())
      create_plots(results_plot())
    })

    # Render the plots
    output$estimatePlot <- shiny::renderPlot({
      shiny::req(plots())
      plots()$estimate_plot
    })

    output$pValuePlot <- shiny::renderPlot({
      shiny::req(plots())
      plots()$p_value_plot
    })


    # When the Submit button is clicked, save the form data
    shiny::observeEvent(input$submit, {
      # Prevent submitting if results are empty
      if(nrow(results_exp()) == 0) {
        shiny::showModal(shiny::modalDialog(
          title = "Error",
          "No results to submit. Please run the simulation first.",
          easyClose = TRUE,
          footer = NULL
        ))
        return()
      }

      # export restuls to database
      #simResults_exp <- appendInputParams_pgsql(results(), input)
      saveData(results())

      # Clear the results after submission
      results_exp(data.frame())

      # Update the responses table with new data
      output$responses <- DT::renderDT({
        loadData()
      }, options = list(pageLength = 5))

    })

    # Conditionally display the Simulation Results header
    output$simulation_results_header <- shiny::renderUI({
      if (nrow(results()) > 0) {
        shiny::h4("Simulation Results")
      } else {
        NULL
      }
    })

    # Download handler for exporting data
    output$downloadBtn <- shiny::downloadHandler(
      filename = function() {
        paste0("Simulation_Results_", Sys.Date(), ".csv")
      },
      content = function(file) {
        # Ensure there is data to download
        #shiny::req(loadData())

        # Write the data to a CSV file
        utils::write.csv(loadData(), file, row.names = FALSE)
      }
    )
  }

  # Run the application
  shiny::shinyApp(ui = ui, server = server)
}
