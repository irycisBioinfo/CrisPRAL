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

ScriptsPATH = paste(getwd(),'/..', '/Scripts', sep = '')
source( str_c( ScriptsPATH, '/Modules/downloadFileModule.R'))
source( str_c( ScriptsPATH, '/Modules/Functions.R' ))

options( shiny.maxRequestSize = 500 * 1024 ^ 2 )

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Simple VCF viewer"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            wellPanel(
                # checkboxInput(inputId = 'multi_sample', label = tags$s('Check for vcf with multiple samples'), value = FALSE),
                
                uiOutput(outputId = 'input_vcf_variable')
            ),  # Well Panel end bracket
            wellPanel(
                
                fileInput(inputId = 'vcf_reference', label = 'Input reference sequence', multiple = FALSE, accept = c('fa','fasta','txt'), placeholder = 'Select file'),
                br(),
                actionButton(inputId = 'run_vcf_reader', label = 'Load')
                
            )
        ),  # SidebarPanel end bracket

        # Show a plot of the generated distribution
        mainPanel(
            h3('Variant Calling File:'),
            DTOutput(outputId = 'vcf_table'),
            h3('Common Operations:'),
            fluidRow(
                column(width = 3,
                       wellPanel(
                           
                           h4('Download selected variants'),
                           checkboxInput(inputId = 'select_all_variants', label = 'Select all', value = FALSE),
                           downloadButtonModule("download_variants", "Download"),
                           actionButton("browser", "browser")
                           
                       ) # Well Panel end bracket
                       ),# Column end bracket
                column(width = 6,
                       wellPanel(
                           
                            h4('Variants'),
                            textOutput(outputId = 'display_variant', inline = TRUE)
                           
                       )) # Well Panel end bracket
                       ),
            fluidRow(
                column(width = 5,
                wellPanel(
                h4('Recalculate Alternate Frequency (AF) values'),
                textOutput(outputId = 'AF_explanation'),
                br(),
                actionButton('Recalculate_AF',label = 'Recalculate AF', icon = icon("redo-alt"))
                
            )))# fluidRow & wellPanel & column end bracket
            

        ) # Main Panel end bracket
    ) # SideBar layout end bracket
) # fluidPage end bracket

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    datos <- reactiveValues()
    
    observeEvent(input$run_vcf_reader, {
        
        vcf_r <- reactive({read.vcfR(input$vcf_file$datapath)})
        datos$vcf_tidy <- vcfR2tidy(vcf_r())
        
        vcf.meta <- reactive({datos$vcf_tidy$meta})
        vcf.table <- reactive({datos$vcf_tidy$fix})
        vcf.gt <- reactive({datos$vcf_tidy$gt})
        
        vcf_reference <- reactive({
            if(is.null(input$vcf_reference)){return('no-reference')
            }else{return(input$vcf_reference$name)}
        })
        
        reference_type <- reactive({return(str_sub(vcf_reference(), 1,2))})
        reference_ID <- reactive({return(str_split(vcf_reference(), '\\.')[[1]][1])})
        
        output$vcf_table <- renderDT(extensions = c('Buttons','FixedColumns','Scroller','ColReorder'), 
                                     options = list(dom = 'Bfrtip', 
                                                    scrollX = TRUE, fixedColumns = list(leftColumns = 2),
                                                    deferRender = TRUE, scrollY = 400, scroller = TRUE,
                                                    colReorder = TRUE,
                                                    buttons = list('copy', list(extend = 'collection',
                                                                               buttons = c('csv', 'excel'),
                                                                               text = 'Download'),
                                                                   'colvis', list(extend = 'colvisGroup',
                                                                                    text = 'INFO',
                                                                                    show = c( 1:6 ),
                                                                                     hide = c( 7:51 )))),
                                     selection = 'multiple',{
                                    vcf.table() %>% select(-ChromKey)
                                    }) #output end brackets
        
        datos$reference_type_ab <- switch(reference_type(), 'NF' = 'c', 'NC' = 'g', 'NP' = 'p', 'NG' = 'g', 'NR' = 'c', 'LR' = 'g')
        
        if(!is.null(datos$reference_type_ab)){
            output$display_variant <- renderText(str_c(reference_ID(),':',
                                                       datos$reference_type_ab,'.',
                                                       vcf.table()[input$vcf_table_row_last_clicked,]$POS,
                                                       vcf.table()[input$vcf_table_row_last_clicked,]$REF,'>',
                                                       vcf.table()[input$vcf_table_row_last_clicked,]$ALT))
        }else{
            output$display_variant <- renderText(str_c(vcf.table()[input$vcf_table_row_last_clicked,]$POS,
                                                       vcf.table()[input$vcf_table_row_last_clicked,]$REF,
                                                       '>',vcf.table()[input$vcf_table_row_last_clicked,]$ALT))
        }
        
    })#observe event end brackets
    
    
    vcf_table_proxy <- DT::dataTableProxy("vcf_table")
    observeEvent(input$select_all_variants, {
        
        if (isTRUE(input$select_all_variants)) {
            DT::selectRows(vcf_table_proxy, input$vcf_table_rows_all)
        } else {
            DT::selectRows(vcf_table_proxy, NULL)
        }
        
    }) # observeEvent end bracket


    
    observeEvent(input$browser,{
        
        browser()
        
    })

    variants_to_download <- reactive({
        
        data <- tibble(str_c(reference_ID(),':',
              datos$reference_type_ab,'.',
              vcf.table()[input$vcf_table_rows_selected,]$POS,
              vcf.table()[input$vcf_table_rows_selected,]$REF,'>',
              vcf.table()[input$vcf_table_rows_selected,]$ALT)) %>% pull()
        
        return(data)
        
    })
    
############# INPUT VCF UI #######################
    #################################.
    
    output$input_vcf_variable <- renderUI({
        # Change multiple = FALSE for multiple = input$multi_sample when ready.
        fileInput(inputId = 'vcf_file', label = 'Input vcf', multiple = FALSE, accept = c('.vcf','.tsv'), placeholder = 'Select vcf file')
        })
    
############# DOWNLOAD FILES #######################
    #################################.
    variant_file_name <- reactive({'variant_list.txt'})
    callModule(downloadFile, 'download_variants', data = variants_to_download, name = variant_file_name)
    
############# recalculate AF #######################
    #################################.
    
    output$AF_explanation <- renderText("Alternate frequency commonly needs to be recalculated, some variant callers that assume diploid or haploid distribution of the data, and in consequence they print 0,0.5 or 1 in the AF column. Despite the real AF value being something in between.")
    observeEvent(input$Recalculate_AF, {datos$vcf_tidy$fix <- datos$vcf_tidy$fix %>% mutate(AF = signif(as.numeric(AO)/DP,digits = 2))})
}



# Run the application 
shinyApp(ui = ui, server = server)
