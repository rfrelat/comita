if (!require(comita)) {
  pak::pak("rfrelat/comita")
  library(comita)
}

# appdir <- here::here("inst/shiny-example/comita/")
# shiny::runApp(appdir)
# # deploy the shinyapp to online server
# rsconnect::deployApp(
#   appDir = appdir,
#   appFiles = rsconnect::listDeploymentFiles(appdir),
#   appName = "Comparative tool for Integrative Trend Analysis",
#   appTitle = "COMITA"
# )
