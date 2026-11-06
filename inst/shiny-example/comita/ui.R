shinyUI(fluidPage(
  titlePanel("Comparative tools for Integrative Trend Analysis"),

  sidebarLayout(
    sidebarPanel(
      conditionalPanel(
        'input.out === "Dataset"',
        fileInput(
          "file1",
          "Choose Data File",
          multiple = FALSE,
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          )
        ),
        radioButtons(
          "sep",
          h6("Separator"),
          choices = c(Tab = "\t", Comma = ",", Semicolon = ";"),
          selected = "\t"
        ),
        hr(),
        selectInput(
          "reg",
          h5("Or select a region:"),
          choices = list("Baltic Sea" = 1, "North Sea" = 2),
          selected = NULL
        ),
        hr(),
        actionButton("choice", "Load dataset")
      ),

      conditionalPanel(
        'input.out === "Results" || input.out === "Comparison"',
        selectInput("met", h5("Method:"), choices = met_choices),
        sliderInput("npc", h5("Number of PC"), min = 2, max = 5, value = 2),
        checkboxInput("log", "Log transform", FALSE)
      ),

      conditionalPanel(
        'input.out === "Results"',
        selectInput(
          "pal",
          h5("Color palette"),
          choices = list("Trafficlight" = 1, "BrBG" = 2, "RdYiBu" = 3)
        ),
        hr()
      ),

      conditionalPanel(
        'input.out === "Results" && input.met === "pca"',
        h5("Horn's parallel analysis"),
        checkboxInput("ran", "Compare with random", FALSE)
      ),

      conditionalPanel(
        'input.out === "Results" && input.met === "pca" && input.ran',
        selectInput(
          "metrand",
          h5("Random model:"),
          choices = list(
            "Random" = "rand",
            "Random + trend" = "trend",
            "Phase structure" = "phase"
          )
        ),
        numericInput("nr", "Repetitions", 100)
      ),

      conditionalPanel(
        'input.out === "Comparison"',
        checkboxGroupInput(
          "multimet",
          h5("Alternative methods"),
          choices = met_choices,
          selected = "mafa"
        ),
        hr(),
        selectInput(
          "npc2",
          h5("Show PC:"),
          choices = list("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5)
        )
      ),

      conditionalPanel(
        'input.out === "Documentation"',
        HTML(
          "<b>Details:</b><br><br>R package comita, v 0.2<br>Last update: 2026-11-06"
        ),
        HTML(
          "<br><br>If you have any question<br>or detect any bug,<br>please contact Romain Frelat:<br>rom.frelat at gmail.com."
        )
      )
    ),

    mainPanel(
      tabsetPanel(
        id = "out",
        tabPanel(
          "Dataset",
          checkboxGroupInput(
            "show_vars",
            label = h5("Select variables"),
            selected = NULL
          )
        ),
        tabPanel(
          "Results",
          plotOutput("distPlot", width = "600px", height = "500px")
        ),
        tabPanel(
          "Comparison",
          plotOutput("compPlot", width = "600px", height = "500px")
        ),
        tabPanel("Documentation", htmlOutput("renderedReport"))
      )
    )
  ),

  hr(),
  HTML(
    "Warning: this is a beta version, use at your own risk.<br>When the screen freezes, please don't make any changes and wait patiently. Computation can take a long time."
  )
))
