library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Comparative tools for Integrative Trend Analysis"),
  #print(h6("v0.4 beta")),
  # Sidebar with a slider input for the number of bins
  sidebarLayout(
   sidebarPanel(
     conditionalPanel(
      'input.out === "Dataset"',
      # Input: Select a file ----
      fileInput("file1", "Choose Data File",
               multiple = FALSE,
               accept = c("text/csv",
                          "text/comma-separated-values,text/plain",
                          ".csv")),
      # Input: Select separator ----
      radioButtons("sep", h6("Separator"),
                  choices = c(Tab = "\t", 
                              Comma = ",",
                              Semicolon = ";"),
                  selected = "\t"),
      hr(),
      selectInput("reg", h5("Or select a region:"),
                  choices = list("Baltic Sea"=1,
                                 "North Sea"=2),
                    # list("Baltic Sea"=1,
                    #              "southern North Sea"=2, 
                    #              "northern North Sea"=3, 
                    #              "English channel"=4,
                    #              "Bay of Biscay"=5), 
                  selected = NULL), 
      hr(),
      actionButton("choice", "Load dataset")
      ),
     
     conditionalPanel(
       'input.out === "Results" | input.out === "Comparison"',
       selectInput("met", h5("Method:"),
                   choices = list("PCA"=1, 
                                  "DFA (very long)"=2, 
                                  "DPCA (long)"=3,
                                  "MAFA"=4, 
                                  "ICA"=5, 
                                  "TSFA"=7, 
                                  "FCA (very long)"=8, 
                                  "LLE"=9)),
       sliderInput("npc", h5("Number of PC"),
                   min = 2, max = 5, value = 2),
       checkboxInput('log', 'Log transform', FALSE)
     ),
       
     conditionalPanel(
       'input.out === "Results"',
       selectInput("pal", h5("Color palette"),
                   choices = list("Trafficlight"=1,"BrBG"=2, "RdYiBu"=3)),
       hr()
     ),
     conditionalPanel(
       'input.out === "Results" & input.met === "1"',
       checkboxInput('ran', 'Compare with random', FALSE)
     ), 
     conditionalPanel(
       'input.out === "Results" & input.met === "1" & input.ran',
       checkboxInput('ar', 'Similar AR', FALSE),
       numericInput("nr", "Repetitions", 100)
     ),
     
     conditionalPanel(
       'input.out === "Comparison"',
       checkboxGroupInput("multimet", 
                          h5("alternative methods"), 
                          choices = list("PCA"=1,
                                         "DFA (very long)"=2, 
                                         "DPCA (long)"=3,
                                         "MAFA"=4, 
                                         "ICA"=5, 
                                         "TSFA"=7, 
                                         "FCA (very long)"=8, 
                                         "LLE"=9),
                          selected = c(7)),
       hr(),
       selectInput("npc2", h5("Show PC:"),
                   choices = list("1"=1, 
                                  "2"=2, 
                                  "3"=3,
                                  "4"=4, 
                                  "5"=5))
       ),
       conditionalPanel(
       'input.out === "Documentation" ',
       HTML("<b>Details: </b> <br> <br> R package comita, v 0.1 <br> Last update : 24.05.2019"),
       HTML("<br><br> If you have any question <br> or detect any bug, <br> please contact R. Frelat: <br> rfrelat(at)yahoo.com.")
       )
    ),
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        id = 'out',
        tabPanel("Dataset", 
                 checkboxGroupInput("show_vars", label = h5("Select variables"),
                                    selected = NULL),
                 verbatimTextOutput("loadread")
                 ),
        tabPanel("Results", plotOutput("distPlot", width="600px", height = "500px")),
        tabPanel("Comparison", plotOutput("compPlot", width="600px", height = "500px")),
        tabPanel("Documentation", htmlOutput("renderedReport"))
      )
    )
  ),
  hr(),
  HTML("Warning: this is a beta-version, use at your own risk."),
  HTML("<br> When screen freeze, please don't make any change and wait patiently. Computation can take a long time.")
))
