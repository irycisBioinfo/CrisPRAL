library(shiny)
library(DT)
require(xlsx)
require(tidyverse)


#ToDo:
#Download the filtered file!

source('Modules/downloadCSVModule.R')

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Lazy Panel Filter"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      navbarPage(title = p(strong('Upload Files')),
        tabPanel(p('Unfiltered File'),
                 br(),
                 # Input: Select a file ----
                 fileInput("file1", "Choose CSV File",
                           multiple = TRUE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv")),
                 
                 # Horizontal line ----
                 tags$hr(),
                 
                 # Input: Checkbox if file has header ----
                 checkboxInput("header", "Header", TRUE),
                 
                 # Input: Select separator ----
                 radioButtons("sep", "Separator",
                              choices = c(Comma = ",",
                                          Semicolon = ";",
                                          Tab = "\t"),
                              selected = ","),
                 
                 # Input: Select quotes ----
                 radioButtons("quote", "Quote",
                              choices = c(None = "",
                                          "Double Quote" = '"',
                                          "Single Quote" = "'"),
                              selected = '"'),
                 
                 # Horizontal line ----
                 tags$hr(),
                 
                 # Input: Select number of rows to display ----
                 radioButtons("disp", "Display",
                              choices = c(Head = "head",
                                          All = "all"),
                              selected = "all")
                 
                 ),
        tabPanel(p('Panel'),
                 br(),
                 # Input: Select a file ----
                 fileInput("file2", "Choose CSV File",
                           multiple = TRUE,
                           accept = c("text/csv",
                                      "text/comma-separated-values,text/plain",
                                      ".csv",
                                      ".xlsx")),
                 
                 # Horizontal line ----
                 tags$hr(),
                 
                 # Input: Checkbox if file has header ----
                 checkboxInput("header2", "Header", FALSE),
                 
                 # Input: Select separator ----
                 radioButtons("sep2", "Separator",
                              choices = c(Comma = ",",
                                          Semicolon = ";",
                                          Tab = "\t"),
                              selected = ","),
                 
                 # Input: Select quotes ----
                 radioButtons("quote2", "Quote",
                              choices = c(None = "",
                                          "Double Quote" = '"',
                                          "Single Quote" = "'"),
                              selected = '"'),
                 
                 # Horizontal line ----
                 tags$hr(),
                 
                 # Input: Select number of rows to display ----
                 radioButtons("disp2", "Display",
                              choices = c(Head = "head",
                                          All = "all"),
                              selected = "all")
                 
        ))
      ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      navbarPage(title = p(strong('Display:')), id = 'Display',
        tabPanel(p('Unfiltered File'),
          # Output: Data file ----
          DTOutput("Unf_file")
          ),
        tabPanel(p('Panel'),
          DTOutput('Genes')
                 ),
        tabPanel(p(strong('Filtered File')), value = 'Filtered_File',
                 fluidRow(column(3, actionButton('Load_filtered_file', 'Start'), uiOutput('react_download_button'))),
                 br(),
                 DTOutput('Fil_file'),
                 br(),
                 fluidRow(
                 column(3, uiOutput('Manual_Options')),
                 column(3, uiOutput('chooseOrg'), uiOutput('chooseDB'), uiOutput('Link'))))
      
    ))
)
)




# Define server logic to read selected file ----
server <- function(input, output) {
  URLdata <- reactiveValues()
  DataTable <- reactiveValues()
  
  output$Unf_file <- renderDT({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    DataTable$UnfilteredDF <- read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)
    
    if(input$disp == "head") {
      return(head(DataTable$UnfilteredDF))
    }
    else {
      return(DataTable$UnfilteredDF)
    }

    
  })
  
  output$Genes <- renderDT({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file2)
    
    # system(paste('./Lazy_Panel_Filter.sh',input$file1$name,input$file2$name,sep=' '))
    
    if (substring(text = input$file2$name, nchar(input$file2$name)-4,nchar(input$file2$name))=='.xlsx'){
      DataTable$PanelDF <- read.xlsx(input$file2$datapath,1,
                     header = input$header2)}
    
    if(substring(text = input$file2$name, nchar(input$file2$name)-3,nchar(input$file2$name))=='.csv'){
      DataTable$PanelDF <- read.csv(input$file2$datapath,
                     header = input$header2,
                     sep = input$sep2,
                     quote = input$quote2)}
                     
    
    if(input$disp2 == "head") {
      return(head(DataTable$PanelDF))
    }
    else {
      return(DataTable$PanelDF)
    }
  })
  
  #filename <- reactive({paste(substring(text = input$file1$name, first = 1, last = nchar(input$file1$name)-3), 'filtered', '.csv', sep='')})
  filename <- reactive({paste(substring(text = input$file1$name, first = 1, last = nchar(input$file1$name)-3), 'filtered', sep='')})
 
  observeEvent(input$Load_filtered_file, {

    for (i in 1:length(unlist(DataTable$PanelDF))){
      
      dt_temp = filter(DataTable$UnfilteredDF, str_detect(DataTable$UnfilteredDF[,'RegionID'],as.character(DataTable$PanelDF[i,])))
      
      if (i == 1){ dt_filtered = dt_temp
      
      }else{
        
        dt_filtered = bind_rows(dt_filtered, dt_temp)
        
      }
    }
    DataTable$FilteredDT <- dt_filtered
    
    output$Fil_file = renderDT({
      
      if(input$disp == "head") {
        return(head(DataTable$FilteredDT))
      }
      else {
        return(DataTable$FilteredDT)
      }}
    
   ,selection = 'single')
    
  })
  
  output$react_download_button <- renderUI({ req(input$Load_filtered_file)
                                              br()
                                              downloadButtonModule(id = 'downloadFile1', 'Download CSV')})
  
  callModule(downloadCSV, id = 'downloadFile1', DataTable$FilteredDT, filename())
  
  observeEvent(input$Fil_file_rows_selected, {
    URLdata$chr <- paste('chr=',substring(DataTable$FilteredDT[input$Fil_file_rows_selected,1],4,nchar(as.character(DataTable$FilteredDT[input$Fil_file_rows_selected,1]))),sep = '')
    URLdata$from <- paste('from=',DataTable$FilteredDT[input$Fil_file_rows_selected,2],sep='')
    URLdata$to <- paste('to=',DataTable$FilteredDT[input$Fil_file_rows_selected,3],sep='')
  })
    
  output$chooseOrg <- renderUI({req(input$Fil_file_rows_selected)
    
  wellPanel(selectInput('Organism',
                                    label = 'Select Organism:',
                                    choices = c('Human', 'Mus musculus')))})
  
  output$chooseDB <- renderUI({req(input$Organism)
                      wellPanel(
                        if (input$Organism == 'Human'){
                          selectInput('DB',
                                      label = 'Select Database:',
                                      choices = c('GRCH38.p12' = 'GCF_000001405.38',
                                                  'GRCH37.p13' = 'GCF_000001405.25',
                                                  'HuRef' = 'GCF_000002125.1',
                                                  'GMH1_1.1' = 'GCF_000306695.2'))}else{
                          selectInput('DB',
                                      label = 'Select:',
                                      choices = c('GRCm38.p6' = 'GCF_000001635.26',
                                                  'Mm_Celera' = 'GCF_000002165.2'))
                                                    
                                                    })})
  
  output$Link <- renderUI({req(input$Fil_file_rows_selected)
    
  wellPanel(helpText(a('Genome Data Viewer',
                      href = paste('https://www.ncbi.nlm.nih.gov/genome/gdv/?context=genome',
                                   paste('acc=', input$DB, sep = ''),
                                   URLdata$chr,
                                   URLdata$from,
                                   URLdata$to,
                                   sep='&'), 
                      target="_blank")))})

  output$Manual_Options <- renderUI({req(input$Fil_file_rows_selected)
      
      wellPanel(textInput('input_chr',
                          'Chromosome:',
                          value = substring(DataTable$FilteredDT[input$Fil_file_rows_selected,1],
                                            4,nchar(as.character(DataTable$FilteredDT[input$Fil_file_rows_selected,1])))),
                          wellPanel(numericInput('input_from',
                                      'Gap Start:',
                                      min = 1,
                                      value = DataTable$FilteredDT[input$Fil_file_rows_selected,2])),
                          wellPanel(numericInput('input_to',
                                      'Gap End:',
                                      min = 1,
                                      value = DataTable$FilteredDT[input$Fil_file_rows_selected,3]))
                                        )})
    

  
}
# Run the app ----
shinyApp(ui, server)