#-Executing piece for All_for_One app.


library(shiny)
source('Scripts/global.R')
source('Scripts/uiMain.R', local = TRUE)
source('Scripts/serverMain.R', local = TRUE)

shinyApp(
  ui = uiMain,
  server = serverMain
)