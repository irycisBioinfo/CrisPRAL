#-Executing piece for All_for_One app.


library(shiny)
source('Scripts/global.R')
source('uiMain.R', local = TRUE)
source('serverMain.R', local = TRUE)

shinyApp(
  ui = uiMain,
  server = serverMain
)
