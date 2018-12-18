#Lazy Panel filter ui

library(shiny)

# Define UI for application that draws a histogram
# shinyUI(
  fluidPage(
  
  # App title ----
titlePanel("Lazy Panel Filter"),
        
        # Sidebar layout with input and output definitions ----
  sidebarLayout(
          
          # Sidebar panel for inputs ----
      sidebarPanel(
            navbarPage(title = p(strong('Upload Files')), id = 'SidebarDisp',
               tabPanel(p('Unfiltered File'), value = 'Unf_file',
                        br(),
                        # Input: Select a file ----
                        fileInput("file1", "Choose CSV File",
                            multiple = TRUE,
                            accept = c("text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv")),
                                
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
                                
                                tags$hr(),
                                
                                # Input: Select number of rows to display ----
                                radioButtons("disp", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "all")
                                
                       ),
               tabPanel(p('Panel'),value = 'Panel',
                        br(),
                        # Input: Select a file ----
                        fileInput("file2", "Choose CSV File",
                            multiple = TRUE,
                            accept = c("text/csv",
                                       "text/comma-separated-values,text/plain",
                                       ".csv",
                                       ".xlsx")),
                                
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
        navbarPage(title = p(strong('Display:')), id = 'MainDisp',
              tabPanel(p('Unfiltered File'),value = 'Unf_file',
                      # Output: Data file ----
                      DTOutput("Unf_file")
              ),
              tabPanel(p('Panel'),value = 'Panel',
                      DTOutput('Genes')
              ),
              tabPanel(p(strong('Filtered File')), value = 'Filtered_File',
                      fluidRow(column(3, 
                                  actionButton('Load_filtered_file', 'Start'),
                                  uiOutput('react_download_button'))),
                      br(),
                      DTOutput('Fil_file'),
                      br(),
                      fluidRow(
                        column(3, uiOutput('Manual_Options')),
                        column(3, uiOutput('chooseOrg'), uiOutput('chooseDB'), 
                               uiOutput('Link'))))
             
              ))
        )
)
  # )
