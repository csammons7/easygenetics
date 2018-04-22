#' Launch Application
#'
#' Implementing the function launch_app() allows the user to launch the app.
#' This app includes five functions.
#'
#' @example launch_app()
#' @export

launch_app <- function() {
  #shiny::runApp("easygenetics", "easygenetics")
  #runApp('R')
   shiny::runApp(system.file('shinyApp', package='easygenetics'))
}
