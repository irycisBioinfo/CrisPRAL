#CrisPRAL ui
#---UI----

#shinyUI(
fluidPage(# Application title
 useShinyjs(),
 tags$head(
  tags$style( #-To scale and center progression bar.
   HTML(".shiny-notification {
        height: 100px;
        width: 800px;
        position:fixed;
        top: calc(50% - 50px);;
        left: calc(50% - 400px);;
        }
        " ))
   ),
 
 # tags$style(type = 'text/css', '#Target_Location {background-color: rgba(255,0,0,0.40);}'), #Example of css editing line
 
 titlePanel( title=(div("Mosaic Finder", img(src="dna_free.png")))),
 
 # Sidebar ----
 #with a slider input for number of bins
 sidebarLayout(
  sidebarPanel(
   checkboxInput("crispr-cas", p(strong('Crispr-Cas9 Analysis')), value = TRUE),
   wellPanel(
    checkboxInput("single_end", p(strong("Single end data"))),
    
    h4("Input Files"),
    fileInput("R1", "Reads 1", accept = c('.fastq', '.fastq.gz')),
    conditionalPanel(
     condition = "input.single_end == false",
     fileInput("R2", "Reads 2", accept = c('.fastq', '.fastq.gz'))),
    fileInput("Reference", "Reference", accept = c('.fasta', '.fastq', '.txt'))
   ),
   wellPanel(inputReset_UI('reset_Target' , label = 'Target')
   ),
   wellPanel(
    h4("Filter Reads Options"),
    numericInput(
     "MinLength",
     "Minimum reads length",
     value = 100,
     min = 31,
     max = 400,
     step = 10),
    wellPanel(
     #Adapter trimming-----
     navbarPage(title = p(strong('Adapter trimming:')), selected = 'recommended',
                tabPanel(p(strong('None')), value = 'None'),
                tabPanel(p(strong('by length')), value = 'length',
                         numericInput('trimA1', 
                                      "R1 sequence",
                                      min = 0,
                                      value = 0,
                                      max = 125),
                         conditionalPanel(
                          condition = "input.single_end == false",
                          numericInput('trimA2', 
                                       "R2 sequence",
                                       min = 0,
                                       value = 0,
                                       max = 125))),
                tabPanel(p(strong('by sequence (recomended)')), 
                         value = 'recommended',
                         tabsetPanel(id = 'input_type_A',
                                     tabPanel('Text input', value = 'Direct_input_A',
                                              br(),
                                              textInput("R1_Adapter","R1 sequence"),
                                              conditionalPanel(
                                               condition = "input.single_end == false",
                                               textInput("R2_Adapter","R2 sequence")),
                                              br(),
                                              fluidRow(#column(4, actionButton('resetP', 'Reset Input')), #No need for this anymore
                                               column(4, actionButton('uploadAdapters', 'Upload Adapters'))),
                                              br(),
                                              br()),
                                     tabPanel('File input', value = 'File_input_A',
                                              hr(),
                                              fileInput("AdapterFile", 'From file (Fasta format)', accept = '.fasta'),
                                              fluidRow(column(4, actionButton('resetA', 'Reset Input'))),
                                              br(),
                                              br()
                                     )),
                         
                         verbatimTextOutput(outputId = 'checkR1A'),
                         conditionalPanel(
                          condition = "input.single_end == false",
                          verbatimTextOutput(outputId = 'checkR2A')))
                
     )),
    wellPanel(
     #Primer Trimming-----
     navbarPage(title = p(strong('Primer trimming:')), selected = 'recommended',
                tabPanel(p(strong('None')), value = 'None'),
                tabPanel(p(strong("by length")), value = 'length',
                         numericInput('trimP1', 
                                      "R1 Trimming",
                                      min = 0,
                                      value = 0,
                                      max = 125),
                         conditionalPanel(
                          condition = "input.single_end == false",
                          numericInput('trimP2', 
                                       "R2 Trimming",
                                       min = 0,
                                       value = 0,
                                       max = 125))),
                tabPanel(title = p(strong("by sequence (recomended)")), value = 'recommended',
                         tabsetPanel(id = 'input_type_P',
                                     tabPanel("Text input",
                                              value = 'Direct_input_P',
                                              hr(),
                                              textInput('F_Primer1',"5' forward end sequence"),
                                              conditionalPanel(
                                               condition = "input.single_end == false",
                                               textInput('R_Primer2',"5' reverse end sequence")),
                                              br(),
                                              fluidRow(#column(4, actionButton('resetP', 'Reset Input')), #No need for this anymore
                                               column(4, actionButton('uploadPrimers', 'Upload Primers'))),
                                              br(),
                                              br()
                                     ),
                                     tabPanel('File input', value = 'File_input_P',
                                              hr(),
                                              fileInput('PrimerFile', 'From file (Fasta format)', accept = c('.fasta','.txt')),
                                              verbatimTextOutput(outputId = 'Explicative_Text1'),
                                              br(),
                                              br())),
                         h4('Trimmed sequence fragment'),
                         verbatimTextOutput(outputId = 'checkFP1'),
                         verbatimTextOutput(outputId = 'checkFP2'),
                         conditionalPanel(
                          condition = "input.single_end == false",
                          verbatimTextOutput(outputId = 'checkRP1'),
                          verbatimTextOutput(outputId = 'checkRP2'))
                         
                         
                )
     ))),
   #---Well stablished parameters----
   checkboxInput("display_advanced", p(strong("Display advanced options"))),
   conditionalPanel(
    condition = "input.display_advanced == true",
    wellPanel(checkboxInput(inputId = "reverse_complement", value = TRUE,
                            p(strong("Check for reverse complement sequences of the reference"))),
              conditionalPanel(condition = "input.
reverse_complement == true", 
                               numericInput(inputId = 'primer_error_rate', 
                                            value = 0.15,
                                            label = "Primer error tolerance",
                                            min = 0,
                                            max = 1))),
     wellPanel(selectInput("general_allType",
                           "Alignment method",
                           c("global", "overlap"),
                           selected = "overlap",
                           ),
               numericInput(inputId = 'gap_open', value = 10, label = 'Alignment gap open penalty', min = 0),
               numericInput(inputId = 'gap_extend', value = 0.5, label = 'Alignment gap extension penalty', min = 0)),
     conditionalPanel(
      condition = "input.single_end == false",    
      wellPanel(
        h4("FastQ-Join parameters"),
        sliderInput(
         "Np",
         "Maximun difference (%)",
         min = 1,
         max = 100,
         value = 0
        ),
        sliderInput(
         "Nm",
         "Minimun overlap (nucleotides)",
         min = 4,
         max = 300,
         value = 10
        )
     )),
    wellPanel(
     h4("Clustering Parameter"),
     sliderInput(
      "cov",
      "Coverage (default: 1)",
      min = 0 ,
      max = 1,
      value = 1
     ),
     sliderInput(
      "identity",
      "Identity (default: 1)",
      min = 0,
      max = 1,
      value = 1
     )
    )),
   
   actionButton("Accept", "Run"),
   actionButton("Run_Example", "Run Example"),
   
   br()
  ),
  #---UI Main Panel----
  mainPanel(
   #Clustering Visualization----
   tabsetPanel(id = 'Visualization', selected = 'Clusterization_&_Alignment',
               tabPanel("Clusterization & Alignment", value = "Clusterization_&_Alignment",
                        navbarPage(title = NULL,id = 'Alignment', selected = 'Reference',
                                   tabPanel('Reference', value = 'Reference',
                                            h3("Results"),
                                            DTOutput("tablaR")
                                   ),
                                   
                                   tabPanel('Target', value = 'Target',
                                            h3("Results"),
                                            DTOutput("tablaT")
                                            
                                   )),
                        
                        uiOutput('temp_text'),
                        uiOutput('alignment_with_filtering'),
                        verbatimTextOutput("alignment"),
                        verbatimTextOutput("score"),
                        uiOutput("Display_Unaligned"),
                        conditionalPanel(
                         condition = "input.display_unaligned == true",
                         uiOutput('Unaligned_Fasta')),
                        uiOutput('Blast_Search'),
                        uiOutput('downloadPAIR_file'),
                        br(),
                        uiOutput('downloadFASTA_file'),
                        br(),
                        uiOutput('downloadFASTAS_clstr'),
                        br()
                        # fluidRow(
                        #  uiOutput('downloadMSA'),
                        #  br())
                        # uiOutput('warningClustrUI'))
               ),
               #Graphics----
               tabPanel("Graphics", value = 'Graphics',  # Show a plot of the generated distribution
                        br(),
                        
                        wellPanel(
                         actionButton(inputId =  "GG", "Generate Plots")),
                        fluidRow(
                         column( width = 4,
                                 wellPanel(
                                  h5('Show:'),
                                  checkboxInput("Pie", p(strong("Pie Chart")), value = TRUE),
                                  checkboxInput("MutFreq", p(strong("Mutation Frequencies")), value = TRUE),
                                  checkboxInput("Delloc", p(strong("Deletions Location")), value = TRUE),
                                  checkboxInput("Delsizes", p(strong("Deletion Sizes")), value = TRUE),
                                  checkboxInput("Inloc", p(strong("Insertions Location")), value = TRUE),
                                  checkboxInput("Insizes", p(strong("Insertion Sizes")), value = TRUE)
                                 )),
                         column( width = 5,
                                 conditionalPanel(
                                  condition = "input.Pie == true",
                                  plotlyOutput(outputId = 'Pie_summary'))
                         )),
                        conditionalPanel(
                         condition = "input.MutFreq == true",
                         wellPanel( 
                          h4('Mutation Frequencies'),
                          plotlyOutput(outputId = 'Mut_Freq')
                         )),
                        splitLayout(
                         conditionalPanel(
                          condition = "input.Delloc == true",
                          wellPanel(
                           h4('Deletions Location'),
                           plotlyOutput(outputId = 'Del_loc')
                          )),
                         conditionalPanel(
                          condition = "input.Delsizes == true",
                          wellPanel(
                           h4('Deletion sizes'),
                           plotlyOutput(outputId = 'Del_sizes')
                          ))),
                        splitLayout(
                         conditionalPanel(
                          condition = "input.Inloc == true",
                          wellPanel(
                           h4('Insertions Location'),
                           plotlyOutput(outputId = 'In_loc')
                          )),
                         conditionalPanel(
                          condition = "input.Insizes == true",
                          wellPanel(
                           h4('Insertion sizes'),
                           plotlyOutput(outputId = 'In_sizes')
                          ))),
                        
                        
                        column(width = 4,
                               textOutput(p(em(outputId = 'Dump'))),
                               textOutput(outputId = 'Dump')
                        )),
               # Download Report----
               tabPanel('Download', value = 'Download',
                        navbarPage(title = NULL,id = 'Alignment2', selected = 'Reference',
                          tabPanel('Reference', value = 'Reference',
                            br(),
                            wellPanel(
                             # h3("Table"),
                             h4("Select entries to be downloaded. Ten entries will be selected by default."),
                             checkboxInput('all_clusters', 'Select all clusters'),
                             hr(),
                             DTOutput("tablaD"),
                             hr(),
                             h4("Selected rows will appear in the table below:"),
                             br(),
                             DTOutput("PrintD"))
                            ),
                          tabPanel('Target', value = 'Target',
                          br(),
                          wellPanel(
                            # h3("Table"),
                            h4("Select entries to be downloaded. Ten entries will be selected by default."),
                            checkboxInput('all_clusters', 'Select all clusters'),
                            hr(),
                            DTOutput("tablaTD"),
                            hr(),
                            h4("Selected rows will appear in the table below:"),
                            br(),
                            DTOutput("PrintTD"))
                          )),
                        fluidRow(
                          conditionalPanel("output.tablaD",
                                           downloadButton("downloadReport", 
                                                          "Download Report"),
                                           downloadButtonModule("downloadFile2", 
                                                                'Download CSV'),
                                           downloadButtonModule("downloadFASTAS", 
                                                                "Download FASTAS")))
                        )
   )
  )))
# )
