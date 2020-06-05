library(shiny)
library(DT)
require(xlsx)


#ToDo:
#Download the filtered file!

# Define UI for data upload app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Lazy Panel Filter"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      navbarPage(title = p(strong('Upload Files')),
        tabPanel(p(em('Unfiltered File')),
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
        tabPanel(p(em('Panel')),
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
        tabPanel(p(em('Unfiltered File')),
          # Output: Data file ----
          DTOutput("Unf_file")
          ),
        tabPanel(p(em('Panel')),
          DTOutput('Genes')
                 ),
        tabPanel(p(strong('Filtered File')), value = 'Filtered_File',
                 actionButton('Load_filtered_file', 'Start'),
                 br(),
                 DTOutput('Fil_file'),
                 br(),
                 uiOutput('Manual_Options'),
                 uiOutput('GenomeNCBI_Search')

        )
      
    ))
)
)


# Define server logic to read selected file ----
server <- function(input, output) {
  
  DataTables <- reactiveValues()
  
  output$Unf_file <- renderDT({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    df <- read.csv(input$file1$datapath,
                   header = input$header,
                   sep = input$sep,
                   quote = input$quote)
    
    if(input$disp == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
    
  })
  
  output$Genes <- renderDT({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file2)
    
    # system(paste('./Lazy_Panel_Filter.sh',input$file1$name,input$file2$name,sep=' '))
    
    if (substring(text = input$file2$name, nchar(input$file2$name)-4,nchar(input$file2$name))=='.xlsx'){
      df <- read.xlsx(input$file2$datapath,1,
                     header = input$header2)}
    
    if(substring(text = input$file2$name, nchar(input$file2$name)-3,nchar(input$file2$name))=='.csv'){
      df <- read.csv(input$file2$datapath,
                     header = input$header2,
                     sep = input$sep2,
                     quote = input$quote2)}
                     
    
    if(input$disp2 == "head") {
      return(head(df))
    }
    else {
      return(df)
    }
  })
  
  filename <- reactive({paste(substring(text = input$file1$name, first = 1, last = nchar(input$file1$name)-3), 'filtered', '.csv', sep='')})
 
  observeEvent(input$Load_filtered_file, {
    
    system(paste('./Lazy_Panel_Filter.sh',input$file1$name,input$file2$name,sep=' '))
    
    DataTables$FilteredDT <- read.csv(paste('./',filename(),sep=''),
                                      header = input$header,
                                      sep = input$sep,
                                      quote = input$quote)
    
    output$Fil_file = renderDT({
      
      if(input$disp == "head") {
        return(head(DataTables$FilteredDT))
      }
      else {
        return(DataTables$FilteredDT)
      }}
    
   ,selection = 'single')
    
  })
  
  observeEvent(input$Fil_file_rows_selected, {
    chr=paste('chr=',substring(DataTables$FilteredDT[input$Fil_file_rows_selected,1],4,nchar(as.character(DataTables$FilteredDT[input$Fil_file_rows_selected,1]))),sep = '')
    acc=paste('acc=','GCF_000001635.26',sep = '')
    from=paste('from=',DataTables$FilteredDT[input$Fil_file_rows_selected,2],sep='')
    to=paste('to=',DataTables$FilteredDT[input$Fil_file_rows_selected,3],sep='')
    
    output$GenomeNCBI_Search <- renderUI({column(3, wellPanel(tabsetPanel(
                                                                tabPanel('Organism',
                                                                         br(),
                                                                          selectInput('Org',
                                                                          label = 'Select:',
                                                                          choices = c("Mus Musculus" = "GCF_000001635.26",
                                                                            "Human" = "GCF_000001405.38"))),
                                                                tabPanel('Database',
                                                                         br(),
                                                                         selectInput('DB',
                                                                                     label = 'Select:',
                                                                                     choices = c('hg19' = 'hg19',
                                                                                                 'hg38' = 'hg38'))))),
                                                    wellPanel(helpText(a('Genome Data Viewer', 
                                                                    href = paste('https://www.ncbi.nlm.nih.gov/genome/gdv/?context=genome',
                                                                                 acc,
                                                                                 chr,
                                                                                 from,
                                                                                 to,
                                                                                 sep='&'), 
                                                                    target="_blank"))))})
    
    output$Manual_Options <- renderUI({column(3, wellPanel(textInput('input_chr',
                                                                     'Chromosome:',
                                                                     value = substring(DataTables$FilteredDT[input$Fil_file_rows_selected,1],
                                                                                       4,nchar(as.character(DataTables$FilteredDT[input$Fil_file_rows_selected,1])))),
                                                  wellPanel(numericInput('input_from',
                                                                      'Gap Start:',
                                                                      min = 1,
                                                                      value = DataTables$FilteredDT[input$Fil_file_rows_selected,2])),
                                                  wellPanel(numericInput('input_to',
                                                                      'Gap End:',
                                                                      min = 1,
                                                                      value = DataTables$FilteredDT[input$Fil_file_rows_selected,3]))
                                                  ))})
    
  })
  
}
# Run the app ----
shinyApp(ui, server)