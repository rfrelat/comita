library(comita)
require(ForeCA)
require(ade4)
require(MARSS)
require(lle)
require(tsfa)
require(fractal)
require(RColorBrewer)
require(vegan)
require(freqdom)
require(fastICA)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  loaddata <- eventReactive(input$choice, {
    #req(input$file1, input$reg)
    if (is.null(input$file1)){
      if(input$reg ==1){
        tab <- read.csv("baltic_gb.txt", sep="\t", row.names = 1)
      }
      if(input$reg ==2){
        tab <- read.csv("north_sea_offshore.txt", sep="\t", row.names = 1)
      } 
    } else {
      req(input$file1)
      # when reading semicolon separated files,
      tryCatch(
        {
          tab <- read.csv(input$file1$datapath, sep = input$sep,
                          header = TRUE, row.names=1)
        },
        error = function(e) {
          stop(safeError(e))
        })
    }
    
    #vars <- names(tab)
    # Update select input immediately after clicking on the action button. 
    updateCheckboxGroupInput(session, "show_vars",
                             choices=names(tab), selected = names(tab))
    return(tab)
  })
  
  output$loadread <- renderText({ 
    tab <- loaddata()
    tab <- subset(tab, select = input$show_vars)
    return("")
  })
  
  output$distPlot <- renderPlot({
    ##Load the dataset
    tab <- loaddata()
    tab <- subset(tab, select = input$show_vars)
    lab <- short(names(tab))#substr(names(tab), 0, 7)
    yr <- row.names(tab) #years
    tab <- as.matrix(tab)
    
    ##Multivariate method
    npc <- input$npc
    
    metconv <- c("pca", "dfa", "dpca", "mafa", 
                 "tsfa", "fca", "lle", "mds")
    mvar <- ita(tab, npc = input$npc, met = metconv[as.numeric(input$met)], 
                sca = TRUE, logt = input$log)
    
    #Change PC to have negative trends
    mvar <- negPC(mvar)
    
    ##Visualization
    layout(matrix(1:4, ncol=2), width = c(1.5,1))
    pal <- list("1"=c("springgreen4","chartreuse3","yellow","darkgoldenrod1", "red"), 
                "2"=c("#a6611a", "#dfc27d", "#f5f5f5", "#80cdc1", "#018571"),
                "3"=c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6"))
    col <- unlist(pal[input$pal])
    #1. Heatmap
    if (!is.null(mvar$co)){ #!input$met%in%c(99)
      par(mar=c(2,4,0,0.1))
      plot_heatmap(mvar, col, shortname=lab)
      mtext("Heatmap", side = 3, line = 0.5, adj = 0)
    } else {
      plot.new()
    }
    
    #2. PC time series
    par(mar=c(0.5,4,1,0.1), xaxs="i", yaxs="i")
    plot_timeseries(mvar)

    #3. Eigen values
    if (!is.null(mvar$eig)){#!input$met%in%c(99)
      par(mar=c(4,2,2,0.1))
      
      #compare to random
      if(input$ran & input$met==1){
        randeig <- ita.nrand(tab, nrep = input$nr, met = metconv[as.numeric(input$met)], sca = TRUE,
                             logt = input$log, metrand = input$metrand, npc = input$npc)
        plot_nrand(mvar, randeig)
        # lines(barx$barx, apply(randeig,2,q95), col="red", 
        #       type="b", pch=18, xpd=NA)
        # signif <- barx$perc>apply(randeig,2,q95)
        # star <- ifelse(signif, "*", "")
        # text(x= barx$barx, y= barx$perc + 1, star, xpd=NA, col="red")
      } else {
        barx <- plot_eig(mvar)
      }
    } else (
      plot.new()
    )
    
    #4. Corcircle or dotchart
    if (!is.null(mvar$co) & input$npc>1){
      par(mar=c(4,2,2,0.1))
      plot_var(mvar)
    } else {
      plot.new()
    }
  })
  
  output$compPlot <- renderPlot({
    ##Load the dataset
    tab <- loaddata()
    tab <- subset(tab, select = input$show_vars)
    lab <- short(names(tab))#substr(names(tab), 0, 7)
    yr <- row.names(tab) #years
    tab <- as.matrix(tab)
    
    ##Multivariate method
    npc2 <- as.numeric(input$npc2)
    # if(npc2>input$npc){
    #   input$npc <- npc2
    # }
    met <- unique(c(as.numeric(input$met), as.numeric(input$multimet)))
    mvar <- multita(tab, npc = input$npc, met = metconv[as.numeric(met)],
             sca = TRUE, logt = input$log)
    
    ##Visualization
    layout(matrix(1:4, ncol=2, byrow = TRUE), width=c(4,1))
    #col <- c("black", brewer.pal(9, "Set1")[met[-1]])
    col <- c("black", rainbow(11)[met[-1]])#
    #1.  PC time series
    leg1 <- plot_compts(tab, mvar, col, showleg = FALSE, pc=npc2)
    plot.new()
    legend("center", legend = leg1$leg, lty = 1, col=leg1$col, xpd=NA, bty="n")
    #2. Variables score
    leg2 <- plot_compvar(tab, mvar, col = col, showleg = FALSE, shortname = lab, pc=npc2)
    plot.new()
    legend("center", legend = leg2$leg, pch = 16, col=leg2$col, xpd=NA, bty="n")
  })
  
  output$renderedReport <- renderUI({           
    includeMarkdown('Documentation.md')           
  })
})
