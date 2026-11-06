# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  raw_data <- eventReactive(input$choice, {
    if (is.null(input$file1)) {
      file <- switch(
        input$reg,
        "1" = "baltic_gb.txt",
        "2" = "north_sea_offshore.txt"
      )
      validate(need(
        file.exists(file),
        paste("Cannot find", file, "in the app folder")
      ))
      read.csv(file, sep = "\t", row.names = 1)
    } else {
      read.csv(
        input$file1$datapath,
        sep = input$sep,
        header = TRUE,
        row.names = 1
      )
    }
  })

  # Fill the variable list whenever a dataset is loaded
  observeEvent(raw_data(), {
    vars <- names(raw_data())
    updateCheckboxGroupInput(
      session,
      "show_vars",
      choices = vars,
      selected = vars
    )
  })

  dataset_vars <- reactive({
    tab <- raw_data()
    req(input$show_vars)
    validate(need(
      length(input$show_vars) >= 2,
      "Select at least two variables"
    ))
    as.matrix(tab[, input$show_vars, drop = FALSE])
  })

  mvar <- reactive({
    tab <- dataset_vars()
    validate(need(
      ncol(tab) >= input$npc,
      "Need at least as many variables as principal components"
    ))
    withProgress(message = "Running ITA...", value = 0.5, {
      negPC(ita(
        tab,
        npc = input$npc,
        met = input$met,
        sca = TRUE,
        logt = input$log
      ))
    })
  })

  randeig <- reactive({
    req(input$ran, input$met == "pca")
    ita.nrand(
      dataset_vars(),
      nrep = input$nr,
      met = input$met,
      sca = TRUE,
      logt = input$log,
      metrand = input$metrand,
      npc = input$npc
    )
  })

  output$distPlot <- renderPlot({
    tab <- dataset_vars()
    lab <- short(colnames(tab))
    fit <- mvar()

    layout(matrix(1:4, ncol = 2), width = c(1.5, 1))
    pal <- list(
      "1" = c("springgreen4", "chartreuse3", "yellow", "darkgoldenrod1", "red"),
      "2" = c("#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571"),
      "3" = c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6")
    )
    col <- unlist(pal[input$pal])

    # 1. Heatmap
    if (!is.null(fit$co)) {
      par(mar = c(2, 4, 0, 0.1))
      plot_heatmap(fit, col, shortname = lab)
      mtext("Heatmap", side = 3, line = 0.5, adj = 0)
    } else {
      plot.new()
      text(0.5, 0.5, "No loadings returned for this method")
    }

    # 2. PC time series
    par(mar = c(0.5, 4, 1, 0.1), xaxs = "i", yaxs = "i")
    plot_timeseries(fit)

    # 3. Eigenvalues (or Horn's comparison)
    if (!is.null(fit$eig)) {
      par(mar = c(4, 2, 2, 0.1))
      if (isTRUE(input$ran) && input$met == "pca") {
        plot_nrand(fit, randeig())
      } else {
        plot_eig(fit)
      }
    } else {
      plot.new()
      text(0.5, 0.5, "No eigen values returned for this method")
    }

    # 4. Correlation circle
    if (!is.null(fit$co) && input$npc > 1) {
      par(mar = c(4, 2, 2, 0.1))
      plot_var(fit)
    } else {
      plot.new()
      text(0.5, 0.5, "No loadings returned for this method")
    }
  })

  output$compPlot <- renderPlot({
    ##Load the dataset
    tab <- dataset_vars()
    lab <- short(colnames(tab))

    ##Multivariate method
    npc2 <- as.numeric(input$npc2)

    met <- unique(c(input$met, input$multimet))
    mvar <- multita(
      tab,
      npc = input$npc,
      met = met,
      sca = TRUE,
      logt = input$log
    )

    ##Visualization
    layout(matrix(1:4, ncol = 2, byrow = TRUE), width = c(4, 1))
    #col <- c("black", brewer.pal(9, "Set1")[met[-1]])
    col <- c("black", rainbow(11)[met[-1]]) #
    #1.  PC time series
    leg1 <- plot_compts(tab, mvar, col, showleg = FALSE, pc = npc2)
    plot.new()
    legend(
      "center",
      legend = leg1$leg,
      lty = 1,
      col = leg1$col,
      xpd = NA,
      bty = "n"
    )
    #2. Variables score
    leg2 <- plot_compvar(
      tab,
      mvar,
      col = col,
      showleg = FALSE,
      shortname = lab,
      pc = npc2
    )
    plot.new()
    legend(
      "center",
      legend = leg2$leg,
      pch = 16,
      col = leg2$col,
      xpd = NA,
      bty = "n"
    )
  })

  output$renderedReport <- renderUI({
    includeMarkdown('Documentation.md')
  })
})
