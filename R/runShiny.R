# runShiny ---------------------------------------
#' Run shiny app made using comita package functions
#' 
#' @export
runShiny <- function() {
  appDir <- system.file("shiny-example", "comita", package = "comita")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `comita`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}