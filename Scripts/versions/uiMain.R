#All For One MAIN ui, RunApp from here

library(shiny)

uiMain <- fluidPage(
 
  # Application Theme
 
 theme = shinytheme("simplex"), 
  useShinyjs(), #To allow JavaScript customization using RCode
# Application title
  titlePanel("All In One"),

  navbarPage( 'Select App', id = 'app',
             tabPanel( 'Mosaic Finder', id = 'Mosaic_Finder',
                      source('./CrisPRAL_ui.R', local = TRUE)$value ),
             tabPanel( 'Lazy Panel Filter', id = 'Lazy_Panel',
                      source('./LazyP_ui.R', local = TRUE)$value ),
             tabPanel('Fastq Quality Analysis', id = 'FastqcRShiny',
                      source('./Fastqc_ui.R', local = TRUE)$value )
  )
)