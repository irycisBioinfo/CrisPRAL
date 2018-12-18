#All For One MAIN ui, RunApp from here

library(shiny)

uiMain <- fluidPage(
  
  # Application title
  titlePanel("All For One"),
  
  navbarPage( 'Select App', id = 'app',
             tabPanel( 'Mosaic Finder', id = 'Mosaic_Finder',
                      source('./CrisPRAL_ui.R', local = TRUE)$value ),
             tabPanel( 'Lazy Panel Filter', id = 'Lazy_Panel',
                      source('./LazyP_ui.R', local = TRUE)$value )
  )
)