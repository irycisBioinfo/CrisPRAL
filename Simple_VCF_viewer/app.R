#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(vcfR)
library(tidyverse)
library(DT)

options( shiny.maxRequestSize = 500 * 1024 ^ 2 )

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            fileInput(inputId = 'vcf_file', label = 'Input vcf', multiple = FALSE, accept = c('.vcf','.tsv'), placeholder = 'Select vcf file'),
            actionButton(inputId = 'run_vcf_reader', label = 'Load')
        ),

        # Show a plot of the generated distribution
        mainPanel(
            DTOutput(outputId = 'vcf_table')
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    vcf_r <- reactive({read.vcfR(input$vcf_file$datapath)})
    vcf_tidy <- reactive({vcfR2tidy(vcf_r())})
    
    vcf.meta <- reactive({vcf_tidy()$meta})
    vcf.table <- reactive({vcf_tidy()$fix})
    vcf.gt <- reactive({vcf_tidy()$gt})
    
    observeEvent(input$run_vcf_reader, {
        
        browser()
        output$vcf_table <- renderDT({
            vcf.table()
        }) #output end brackets

    })#observe event end brackets
}

# Run the application 
shinyApp(ui = ui, server = server)
