#ui_Mosaic_batch

# This script includes the user-interface definition of the app.

###############################################.
## Header ---- 
###############################################.

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
  
  ###############################################.
  ## Mosaic Finder batch ----
  ###############################################.
  
  titlePanel( title=(div("Mosaic Finder BATCH", img(src="dna_free.png")))), # titlePanel end bracket
    mainPanel(width = 9, style="margin-left:4%; margin-right:4%",
      wellPanel(
        # checkboxInput("single_end_batch", p(strong("Single end data"))), # Choose between single end or paired end inputs
        
        h4("Mandatory input Files"),
        wellPanel(
        br(),
        fileInput("dir", "Select .fastq data", accept = c('.zip','.gz','.fastq'), multiple = TRUE),
        br(),br(),
        fileInput("Reference_batch", "Reference", accept = c('.fasta', '.fastq', '.txt'))
        ), # wellPanel end bracket
        h4("Optional files"),
        wellPanel(
          
          ##############################################.
          # Adapter trimming-----
          ##############################################.
        fluidRow(column(5,
          navbarPage(title = p(strong('Adapter trimming:')), selected = 'recommended_batch', id = 'Adapter_trimming_batch',
           tabPanel(p(strong('None')), value = 'None_batch'),
           # tabPanel(p(strong('by length')), value = 'length_batch',
           #  numericInput('trimA1_batch', 
           #               "R1 sequence",
           #               min = 0,
           #               value = 0,
           #               max = 125),
           #  conditionalPanel(
           #    condition = "input.single_end == false",
           #    numericInput('trimA2_batch', 
           #                 "R2 sequence",
           #                 min = 0,
           #                 value = 0,
           #                 max = 125))),
           tabPanel(p(strong('by sequence (recomended)')), 
                    value = 'recommended_batch',
                    tabsetPanel(id = 'input_type_A_batch',
                      # tabPanel('Text input', value = 'Direct_input_A_batch',
                      #    br(),
                      #    textInput("R1_Adapter_batch","R1 sequence"),
                      #    conditionalPanel(
                      #      condition = "input.single_end_batch == false",
                      #      textInput("R2_Adapter_batch","R2 sequence")),
                      #    br(),
                      #    fluidRow(#column(4, actionButton('resetP', 'Reset Input')), #No need for this anymore
                      #      column(4, actionButton('uploadAdapters_batch', 'Upload Adapters', icon = icon('upload')))),
                      #    br(),
                      #    br()),
                      tabPanel('File input', value = 'File_input_A_batch',
                         hr(),
                         fileInput("AdapterFile_batch", 'From file (Fasta format)', accept = '.fasta'),
                         fluidRow(column(4, actionButton('resetA_batch', 'Reset Input'))),
                         verbatimTextOutput(outputId = 'Explicative_TextA_batch'),
                         br(),
                         br()
                      )),
                    h4('Sequence to be trimmed'),
                    verbatimTextOutput(outputId = 'checkR1A_batch'), # Prints adapters for the user to check
                    conditionalPanel(
                     condition = "input.single_end == false",
                     verbatimTextOutput(outputId = 'checkR2A_batch')))
            ) # column end bracket
          ), # navPar end bracket
        
          ##############################################.
          #Primer Trimming-----
          ##############################################.
          column(5,
            navbarPage(title = p(strong('Primer trimming:')), selected = 'recommended', id = "Primer_trimming_batch",
               tabPanel(p(strong('None')), value = 'None'),
               # tabPanel(p(strong("by length")), value = 'length',
               #      numericInput('trimP1_batch', 
               #                   "R1 Trimming",
               #                   min = 0,
               #                   value = 0,
               #                   max = 125),
               #      conditionalPanel(
               #        condition = "input.single_end == false",
               #        numericInput('trimP2_batch', 
               #                   "R2 Trimming",
               #                   min = 0,
               #                   value = 0,
               #                   max = 125))),
               tabPanel(title = p(strong("by sequence (recomended)")), value = 'recommended',
                    tabsetPanel(id = 'input_type_P_batch',
                  #       tabPanel("Text input", value = 'Direct_input_P_batch',
                  #          hr(),
                  #          textInput('F_Primer1_batch',"5' forward end sequence"),
                  #          conditionalPanel(
                  #            condition = "input.single_end_batch == false",
                  #            textInput('R_Primer2_batch',"5' reverse end sequence")),
                  #          br(),
                  #          fluidRow(#column(4, actionButton('resetP', 'Reset Input')), #No need for this anymore
                  #            column(4, actionButton('uploadPrimers_batch', 'Upload Primers', icon = icon('upload')))),
                  #          br(),
                  #          br()
                  # ),
                        tabPanel('File input', value = 'File_input_P_batch',
                           hr(),
                           fileInput('PrimerFile_batch', 'From file (Fasta format)', accept = c('.fasta','.txt')),
                           verbatimTextOutput(outputId = 'Explicative_TextP_batch'),
                           br(),
                           br())),
                        h4('Sequence to be trimmed'),
                        # Prints primers for the user to check
                        verbatimTextOutput(outputId = 'checkFP1_batch'),
                        verbatimTextOutput(outputId = 'checkFP2_batch'),
                        conditionalPanel(
                          condition = "input.single_end_batch == false",
                          verbatimTextOutput(outputId = 'checkRP1_batch'),
                          verbatimTextOutput(outputId = 'checkRP2_batch'))
               ) #tabPanel (recommended) closing bracket
              ) # NavBar Page Primer trimming closing bracket         
            ) # primer column closing bracket
          ) # fluidRow closing bracket
        ), # WellPanel PRIMER trimming closing bracket
        wellPanel(
        fluidRow(
            column(5, 
               h4('Score threshold: '),h5('clusters with score < threshold will be filtered out onto a separate sheet'),br(),
               checkboxInput('detect_theshold', label = ' Auto detect optimun scoring threshold', value = TRUE),
               conditionalPanel(condition = "input.detect_theshold == false",
               numericInput(inputId = 'scoring_threshold', 
                            value = 0,
                            label = "Select a scoring threshold",
                            min = 0,
                            max = 99999)) # conditional panel end bracket
             ), # column end bracket
            column(5,
               h4('Indel interest range: '), h5('Clusters with indels within interest range will be selected and printed onto a different sheet'),
               checkboxInput('interest_range_select', label = 'Click to select min and max bounds', value = FALSE),
               br(),
               conditionalPanel(condition = "input.interest_range_select == true",
                  numericInput(inputId = 'min_interest_range', 
                               value = 0,
                               label = "Select a lower bound",
                               min = 0,
                               max = 99999),
                  numericInput(inputId = 'max_interest_range', 
                               value = 0,
                               label = "Select an upper bound",
                               min = 0,
                               max = 99999)
                      ),
               checkboxInput(inputId = "print_mismatches", label = 'Append mismatches breakdown to final table', value = FALSE),
               checkboxInput(inputId = "print_sequences", label = 'Append representative sequences to final table (beta)', value = FALSE)
            ) # column end bracket
          ) #fluidRow end bracket
        ), #wellPanel end bracket
        
        #####################################################.
        ############ Advanced options refuse to work ########.
        #####################################################.
        
        # wellPanel(
        #   checkboxInput(inputId = "display_advanced", p(strong("Display advanced options"))),
        #   conditionalPanel(
        #     condition = "input.display_advanced == true",
        #     wellPanel(
        #     verbatimTextOutput(outputId = 'advanced_for_the_time_being'))
        #     ) # Conditionalpanel end bracket
        # ), #wellPanel end bracket
        
        wellPanel(
        actionButton("Accept_batch", "Run", icon = icon('caret-right')))
                ) # wellPanel end bracket
    ) # mainPanel end bracket
) # fluidPage end bracket