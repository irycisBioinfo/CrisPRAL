#
# 
# TO DO
# *Indicar cual es la secuancia mas parecida al target.
# *Crear graficas de: * crispresso
#       - % de mutacion por base
#       - % de indel por base
#       - % cambio por base
#   * Crear un report en PDF (con markdown) con las graficas.
# 
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/

#---Libraries----
library(base)
base_append <- base::append
library(shiny)
library(plotly)
library(tidyverse)

if (require('devtools') == FALSE){
  
  install.packages('devtools')
  library(devtools)
  devtools::install_github('ropensci/plotly')
  library(plotly)
  
}

if (require('kableExtra') == FALSE){
  
  install.packages('kableExtra')
  library(kableExtra)
  
}

if (require('rmarkdown') == FALSE){
  
  install.packages('rmarkdown')
  library(rmarkdown)
  
}
if (require('shinyBS') == FALSE){
  
  install.packages('shinyBS')
  library('shinyBS')
  
}

if (require('webshot') == FALSE){
  
  install.packages('webshot')
  library(webshot)
  webshot::install_phantomjs()
  
}

if (!require("processx")) install.packages("processx")

#To install Biostrings:
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")

library(Biostrings) 
library(DT)

#------Setting up work environment-----

appPATH = getwd()

if (substring(appPATH, nchar(appPATH)-8, nchar(appPATH)) != 'CrisPRAL'){
  
  system('find /home/ -name "CrisPRAL" > PATH.txt')
  
  appPATH = read_file('PATH.txt')
  
  system('rm PATH.txt')
  
  CompletePATH = substring(appPATH, 1, nchar(appPATH)-1) # nchar - 1 to remove "/n" due to .txt parsing
  
  setwd(CompletePATH)
  
}else{CompletePATH = appPATH}

#---Modules_Load----

source(paste(CompletePATH, '/../Modules/Functions.R', sep = ''))
source(paste(CompletePATH, '/../Modules/input-reset_module.R', sep = ''))
source(paste(CompletePATH, '/../Modules/downloadFileModule.R', sep = ''))


#---UI Side Panel---- 

ui <- fluidPage(# Application title
  
  tags$head(
    tags$style( #-To scale and center progression bar.
      HTML(".shiny-notification {
              height: 100px;
              width: 800px;
              position:fixed;
              top: calc(50% - 50px);;
              left: calc(50% - 400px);;
            }
           " ))),
  
  titlePanel( title=(div("MOSAIC Finder", img(src="dna_free.png")))),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        checkboxInput("single_end", p(strong("Single end data"))),
        
        h4("Input Files"),
        fileInput("R1", "Reads 1", accept = c('.fastq', '.fastq.gz')),
        conditionalPanel(
          condition = "input.single_end == false",
          fileInput("R2", "Reads 2", accept = c('.fastq', '.fastq.gz'))),
        fileInput("Reference", "Reference", accept = c('.fasta', '.fastq'))
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
          #Adapter trimming is now a dummy button, it doesnt lead anywhere. 
          #Fix with a similar structure to that of Primer filtering
        navbarPage(title = p(strong('Adapter trimming')),
          tabPanel(p(strong('None'))),
          tabPanel(p(strong('by length')),
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
          tabPanel(p(strong('by sequence')),
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
        navbarPage(title = p(strong('Primer trimming:')),
          tabPanel(p(strong('None'))),
          tabPanel(p(strong("by length")),
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
          tabPanel(p(strong(("by sequence"))),
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
                        fileInput('PrimerFile', 'From file (Fasta format)', accept = '.fasta'),
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
          "Coverage",
          min = 0 ,
          max = 1,
          value = 1
        ),
        sliderInput(
          "identity",
          "Identity",
          min = 0,
          max = 1,
          value = 1
        )
      ),
      
      actionButton("Accept", "Run")
    ),
#---UI Main Panel----
      mainPanel(
    
        tabsetPanel(id = 'Visualization', selected = 'Clusterization_&_Alignment',
          tabPanel("Clusterization & Alignment", value = "Clusterization_&_Alignment",
                   
                   h3("Table"),
                   DTOutput("tabla"),
                   h3("Alignment"),
                   selectInput(
                     "alnType",
                     "Alignment method",
                     c("global", "local", "overlap", "global-local", "local-global")
                   ),
                    uiOutput('is.nullTarget'),
                    verbatimTextOutput('Print2'),
                    verbatimTextOutput("alignment"),
                    verbatimTextOutput("score"),
                    uiOutput('Blast_Search')
                   
                   ),
                   
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
                      column( width = 7,
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
                         condition = "input.Insizes == true",
                           wellPanel(
                               h4('Insertion sizes'),
                               plotlyOutput(outputId = 'In_sizes')
                   )),
                       conditionalPanel(
                         condition = "input.Inloc == true",
                           wellPanel(
                             h4('Insertions Location'),
                             plotlyOutput(outputId = 'In_loc')
                   ))),
                   
                    column(width = 4,
                      textOutput(p(em(outputId = 'Dump'))),
                      textOutput(outputId = 'Dump')
                     )),
          tabPanel('Download', value = 'Download',
                   br(),
                    wellPanel(
                      tabsetPanel(
                      tabPanel('Download CSV',
                               hr(),
                               downloadButtonModule("downloadFile2", 'Download CSV')),
                      tabPanel('Download Full Report',
                               hr(),
                               h3("Table"),
                               DTOutput("tabla2"),
                               br(),
                               hr(),
                               verbatimTextOutput(outputId = 'Print'),
                               downloadButton("downloadReport"))
                   )))
                  )
  )))

#---SERVER---- 
# Define server logic required to draw a histogram
server <- function(input, output) {

  #---Reactive values----  
  options(shiny.maxRequestSize = 100 * 1024 ^ 2)
  datos <- reactiveValues()
  charts <- reactiveValues()
  file_upload <- reactiveValues(upload_Adapter_state = NULL)
  file_upload <- reactiveValues(upload_Primer_state = NULL)
    #-Reactive input to allow for compressed files processing.
  R1_file <- reactive({
    if (str_sub(input$R1$datapath, -3) == '.gz'){
      Unzip_file(CompletePATH, input$R1$datapath)}
    else{
      input$R1$datapath}
    })
  R2_file <- reactive({
    if(!is.null(input$R2)){
      if(str_sub(input$R2$datapath, -3) == '.gz'){
        Unzip_file(CompletePATH, input$R2$datapath)}
      else{
        input$R2$datapath}
  }else{
    return('Empty')
  }})
    #-Allows reference Primer trimming.
  Reference <- reactive({ 
    if(file_upload$upload_Primer_state == 'uploaded'){
      T_Refence_dir = paste(str_sub(input$Reference$datapath, 1,-8),'Filtered_Reference.fasta', sep = '')
      system(paste(CompletePATH, '/bin/cutadapt/cutadapt -n 2 -g ', as.character(PrimerR1_5()), ' -a ', as.character(PrimerR2_5()),
                   ' -o ', str_sub(input$Reference$datapath, 1,-8), 'Filtered_Reference.fasta ', input$Reference$datapath, sep = '' ))
      return(T_Refence_dir)
    }
    else{
      return(input$Reference$datapath)
    }
    })

  #----Filtering Parameters declarations----
  Target_FileInput <- callModule(inputReset_server, 'reset_Target')

  output$is.nullTarget <- renderUI({#Generates reactive UI depending on the presence of a Target or not.
    selectInput(
      "alnTo",
      "Alignment to:",
      if (!is.null(Target_FileInput())){c("Reference","Target")}else{c("Reference")})
  })

  #This following set of observeEvents is to allow user to change his mind on uploaded Adapter files.
  observeEvent(c(input$AdapterFile, input$uploadAdapters), {
    file_upload$upload_Adapter_state <- 'uploaded'
    })
  observeEvent(c(input$resetA, input$input_type_A), {
    file_upload$upload_Adapter_state <- 'reset'
    })

  Adapter_Directinput <- eventReactive(input$uploadAdapters,{ c(input$R1_Adapter,input$R2_Adapter) })
  Adapter_FileInput <- reactive({
    if (is.null(file_upload$upload_Adapter_state)) {
      return(NULL)
    } else if (file_upload$upload_Adapter_state == 'uploaded') {
      
      if (input$input_type_A == 'File_input_A'){return(input$AdapterFile)}
      else{return(Adapter_Directinput())}
    }
    else if (file_upload$upload_Adapter_state == 'reset') {
      return(NULL)
    }
  })
  
  Adapter_sequences <- reactive({
    
    if (is.null(Adapter_FileInput())){
      return("Empty")
    }else if (!is.null(Adapter_FileInput())){
      if(input$input_type_A == 'File_input_A'){return(readDNAStringSet(Adapter_FileInput()$datapath))
      }else{
        return(DNAStringSet(Adapter_FileInput()))}
    }})
  
  AdapterR1 <- reactive({
    if (is.null(Adapter_FileInput())){
      return("Empty")
    }
    return(as.character(Adapter_sequences()[1]))
  })
  AdapterR2 <- reactive({
    if (is.null(Adapter_FileInput())){
      return("Empty")
    }
    return(as.character(Adapter_sequences()[2]))
  })
  
  output$checkR1A <- renderPrint(paste("R1 Adapter: ",AdapterR1(), sep = ''))
  output$checkR2A <- renderPrint(paste("R2 Adapter: ",AdapterR2(), sep = ''))
  
  #Same than the Adapter strings but with the Primers.
  observeEvent(c(input$PrimerFile, input$uploadPrimers), {
    file_upload$upload_Primer_state <- 'uploaded'
  }) 
  observeEvent(c(input$resetP, input$input_type_P), {
    file_upload$upload_Primer_state <- 'reset'
  })
  
  Primer_Directinput <- eventReactive(input$uploadPrimers,{ c(input$F_Primer1,input$R_Primer2) })
  Primer_FileInput <- reactive({
    if (is.null(file_upload$upload_Primer_state)) {
      return(NULL)
    } else if (file_upload$upload_Primer_state == 'uploaded') {
      
      if (input$input_type_P == 'File_input_P'){return(input$PrimerFile)}
      else{return(Primer_Directinput())}
        }
     else if (file_upload$upload_Primer_state == 'reset') {
      return(NULL)
    }
  })
  
  #Primer upload through file:
  #Primer sequences are loaded as reactive functions to prevent residual data lingering.
  Primer5_sequences <- reactive({

    if (is.null(Primer_FileInput())){
      return("Empty")
    }else if (!is.null(Primer_FileInput())){
      if(input$input_type_P == 'File_input_P'){return(readDNAStringSet(Primer_FileInput()$datapath))
      }else{
          return(DNAStringSet(Primer_FileInput()))}
      }})
  Primer3_sequences <- reactive({

    if (is.null(Primer_FileInput())){
      return("Empty")
    }else if (!is.null(Primer_FileInput())){
    return(reverseComplement(Primer5_sequences()))
    }})
  
  PrimerR1_5 <- reactive({
    if (is.null(Primer_FileInput())){
      return("Empty")
    }
    return(as.character(Primer5_sequences()[1]))
  })
  PrimerR1_3 <- reactive({
    if (is.null(Primer_FileInput())){
      return("Empty")
    }
    return(as.character(Primer3_sequences()[2]))
  })
  PrimerR2_5 <- reactive({
    if (is.null(Primer_FileInput())){
      return("Empty")
    }
    return(as.character(Primer5_sequences()[2]))
  })
  PrimerR2_3 <- reactive({
    if (is.null(Primer_FileInput())){
      return("Empty")
    } 
    return(as.character(Primer3_sequences()[1]))
  })
  
  output$Explicative_Text1 <- renderText(paste('The file should contain four lines, two', 
                                               'lines with ">" indicating Primer ID and', 
                                               'another two lines with the sequences', 'themselves. Ex:',
                                               '>Primer1', 'ACGT', '>Primer2', 'GTAC',sep = '\n'))
  
  output$checkFP1 <- renderPrint(paste("Read 1, 5' end: ",PrimerR1_5(), sep = ''))
  output$checkFP2 <- renderPrint(paste("Read 1, 3' end: ",PrimerR1_3(), sep = ''))
  output$checkRP1 <- renderPrint(paste("Read 2, 5' end: ",PrimerR2_5(), sep = ''))
  output$checkRP2 <- renderPrint(paste("Read 2, 3' end: ",PrimerR2_3(), sep = ''))

  #----Processing----
  
  observeEvent(input$Accept, {
    
    req(input$R1, input$Reference)
    
    withProgress(message = "Performing analisis", {
      
      command = paste(
        CompletePATH,
        "/pipeline_v6.pl --r1 ",   ### Automatizar el PATH del script -- Doned
        R1_file(),
        " --r2 ",
        R2_file(),
        " --single_end ",
        input$single_end,
        " --min_len ",
        input$MinLength,
        " --Adapter_R1 ",
        as.character(AdapterR1()),
        " --Adapter_R2 ",
        as.character(AdapterR2()),
        " --trimA1 ",
        input$trimA1,
        " --trimA2 ",
        input$trimA2,
        " --F_Primer1 ",
        as.character(PrimerR1_5()),
        " --F_Primer2 ",
        as.character(PrimerR1_3()),
        " --R_Primer1 ",
        as.character(PrimerR2_5()),
        " --R_Primer2 ",
        as.character(PrimerR2_3()),
        " --trimP1 ",
        input$trimP1,
        " --trimP2 ",
        input$trimP2,
        " --Np ", 
        input$Np,
        " --Nm ",
        input$Nm,
        " --cov ",
        input$cov,
        " --id ",
        input$identity,
        collapse = "",
        sep = ""
      )
      
      system(command)#-Executes perl script
      
      datos$fasta = read_tsv("cluster.tsv", col_names = FALSE)
      datos$cluster = read_tsv("cluster.bak.tsv", col_names = FALSE)
      
      if (isEmpty(colnames(datos$fasta))){#-If dataset is completely filtered due to any restriction app crashes
        stopApp(returnValue = 'Pipe was broken')
        }else{
        colnames(datos$fasta) = c("ID", "Sequence")
        colnames(datos$cluster) = c("ClusterN", "Length", "ID", "Identity")
        } #WE NEED TO BE ABLE TO SUPPLY ERROR INFORMATION ASAP
      
      datos$Tabla_raw = datos$cluster %>% 
        group_by(ClusterN) %>% 
        mutate(Abundance = n()) %>% 
        ungroup() %>% 
        mutate(Freq = 100 *Abundance / n()) %>%
        filter(Identity == "*") %>%
        select(ID, Abundance, Freq) %>% 
        arrange(desc(Abundance))
      
      Total_Abundance = sum(datos$Tabla_raw['Abundance'])
      
      datos$Ref = readDNAStringSet(Reference())
      
      if(!is.null(Target_FileInput())){datos$Target = readDNAStringSet(Target_FileInput()$datapath)}
        
      datos$Sequences = readDNAStringSet("cluster")
      datos$aln = pairwiseAlignment(datos$Sequences, datos$Ref, type = "global", gapOpening=10, gapExtension=5) #maybe local
      tmp = data.frame(
        ID = names(datos$Sequences),
        mismatch = nmismatch(datos$aln),
        length = nchar(datos$aln),
        score = score(datos$aln),
        width = width(pattern(datos$aln)),
        start = start(pattern(datos$aln)),
        end = end(pattern(datos$aln)),
        deletions = nindel(datos$aln)@deletion[, 2],
        insertions = nindel(datos$aln)@insertion[, 2]
      )
      tmp = tmp %>% separate(ID, c("ID","kk"), sep =" ") %>% select(-kk)
      
      system("rm R1_* R2_* join* good* bad* cluster* final*"); #To prevent constant cashing on the same data
      
      datos$Tabla = inner_join(datos$Tabla_raw,tmp) %>% mutate(score = round(score,1), Freq = signif(Freq,2)) %>% 
        arrange(desc(Abundance))

        output$tabla = renderDT(datos$Tabla, selection = "single")
        output$tabla2 = renderDT(datos$Tabla) #Table with multiple selection for PDF formatting
        output$Print = renderPrint(datos$Tabla[sort.default(input$tabla2_rows_selected),])
        output$Print2 = renderPrint(paste("Total amount of Reads: ", Total_Abundance, sep = ''))
        
        #-Pie-Chart-processing----
        Other <- datos$Tabla %>% filter(Freq < 1)
        datos$Pie_data <- select(datos$Tabla, ID, Freq) %>% filter(Freq > 1) %>% 
                          add_row(ID = paste('<1% groups', 
                                             '(', count(Other['Freq']) , ')', 
                                             sep = ''), 
                                  Freq = sum(Other['Freq']))
        
        #-Pie-Chart-plotting online
        colors <- c('rgb(211,94,96)', 'rgb(128,133,133)', 'rgb(144,103,167)', 'rgb(171,104,87)', 'rgb(114,147,203)')
          
          PC_online <- plot_ly(datos$Pie_data, labels = ~ID, values = ~Freq, type = 'pie',
                        hoverinfo = 'text',
                        marker = list(colors = colors, line = list(color = '#FFFFFF', width = 1)), 
                        text = ~ID ,showlegend = TRUE) %>%
          layout(title = 'Main Groups',
                 xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                 yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
          
          output$Pie_summary <- renderPrint({print(PC_online)})
          
        #-Pie-Chart for plotting in PDF  
          PC_pdf <- plot_ly(datos$Pie_data, labels = ~ID, values = ~Freq, type = 'pie',
                            marker = list(colors = colors, line = list(color = '#FFFFFF', width = 1)), 
                            showlegend = TRUE) %>%
            layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                   yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))
          
          charts$tmpFilePie <- tempfile(fileext = '.png')
          export(PC_pdf, charts$tmpFilePie)
        
        output$Pie_summary <- renderPlotly({print(PC_online)})
        #----
    })
  
  })
  
  observeEvent(c(input$tabla_rows_selected, input$alnType, input$alnTo), {
    withProgress(message = "Performing analisis", {
      
      req(datos$Tabla)
      req(input$tabla_rows_selected)
      
      S = datos$Tabla %>% filter(ID == datos$Tabla$ID[input$tabla_rows_selected]) %>% left_join(datos$fasta)
      Seq = DNAStringSet(S$Sequence)
      if(input$alnTo == "Reference")
      {
        aln = pairwiseAlignment(datos$Ref, Seq, type = input$alnType)
      }else{
        aln = pairwiseAlignment(datos$Target, Seq, type = input$alnType)
      }
      
      
      comp = character(length(aln@pattern))
      ref = unlist(strsplit(as.character(aln@pattern), split = ""))
      query = unlist(strsplit(as.character(aln@subject), split = ""))
      
      for (i in 1:length(ref))
      {
        if (ref[i] != query[i])
        {
          comp[i] = query[i]
        } else {
          comp[i] = "."
        }
      }
      
      pos <- Position_vector(as.character(aln@pattern), ref) #Generates suitable positioning vector
      
      if (!is.null(input$tabla_rows_selected))
      {
        datos$text = as.data.frame(c(
          paste(pos, sep = '', collapse = ''), #Not working well enough, fails in presence of inserts
          paste(ref, sep = "", collapse = ""),
          paste(query, sep = "", collapse = ""),
          paste(comp, sep = "", collapse = "")
        ))
        
        if(input$alnTo == "Reference")
        {
          rownames(datos$text) = c("Position", "Reference", datos$Tabla$ID[input$tabla_rows_selected], "Comparisson")
        }else{
          rownames(datos$text) = c("Position", "Target", datos$Tabla$ID[input$tabla_rows_selected], "Comparisson")
        }
        
        colnames(datos$text) = c("Alignment")
        output$alignment = renderPrint(datos$text)
        output$score = renderPrint(
          paste(
            "Alignment Score:",
            as.numeric(aln@score),
            "Alignment Length",
            length(ref),
            sep = "  ",
            collapse = ""
          )
        )
      }
      output$Blast_Search <- renderUI({column(3, wellPanel(helpText(a(div(img(src="BLASTn.png")), 
                                                                      href = paste('https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch',
                                                                                   '&JOB_TITLE=Cluster%20',
                                                                                   input$tabla_rows_selected,
                                                                                   '&QUERY=>',
                                                                                   datos$Tabla$ID[input$tabla_rows_selected],
                                                                                   '%0A',
                                                                                   Seq, sep=''), 
                                                                      target="_blank"))))})
    })
  })
  
  observeEvent(c(input$GG), { #---Plotting Graphs

        req(datos$Tabla) # Prevents APP from crashing

        #-----------Setting up the Data----------------
    withProgress(message = "Plotting...", { 
          #Extracting Sequences, identificators, sizes
        VCAbundance = select(datos$Tabla, Abundance)
        VCSeq = list(as.character(datos$aln@pattern))
        VCIDs = select(datos$Tabla, ID)
        SizeList = list(nchar(datos$aln@pattern))
        
        Graph_Data = bind_cols(VCIDs, SizeList, VCSeq, VCAbundance) #Bind Identificator with sequence and size as to not lose them 
                                                                    # when we filter by size
        colnames(Graph_Data) = c('ID', 'Size', 'Sequence', 'Abundance') #Name the tables
        
        Filtered_Graph_Data = filter(Graph_Data, Size < 250*(2-input$cov))#Leave enough margin in order to correct for coverage lax
        Dump = filter(Graph_Data, Size > 250*(2-input$cov))
        
        TableSize = max(Filtered_Graph_Data['Size'])
  
        #----------------Processing-------------------------
          #-Mutation Frequency
        IDsSize = (count(Filtered_Graph_Data))*TableSize
        
        LetterPerPosition = data.frame(ID = rep(Filtered_Graph_Data$ID,TableSize), 
                                       Position = sort(rep(1:TableSize,count(Filtered_Graph_Data))),
                                       Chr = substring(Filtered_Graph_Data$Sequence, 
                                                       sort(rep(1:TableSize,count(Filtered_Graph_Data))), 
                                                       sort(rep(1:TableSize,count(Filtered_Graph_Data)))),
                                       Abundance = rep(Filtered_Graph_Data$Abundance, TableSize))
        
        LetterPerPosition_Total = LetterPerPosition %>% group_by(Position, Chr, Abundance) %>% 
                                                        count(Chr) %>% filter(Chr != '') %>% rename(n = 'Count') %>% 
                                                        transmute(Total = Abundance*Count) %>% ungroup() %>%
                                                        select(Total, Chr, Position) %>% group_by(Position, Total) %>% distinct(Chr)
                                      
        
        LetterPerPosition_Norm = LetterPerPosition_Total %>% ungroup() %>% group_by(Position) %>% 
                                                             mutate_at(vars(Total),funs(./sum(Total))) %>% 
                                                             ungroup(Position) %>% group_by(Chr, Position) %>%
                                                             mutate_at(vars(Total), funs(sum(Total))) %>% 
                                                             group_by(Total, Position) %>% distinct(Chr)
        
          #-Isolating deletions
        Deletions_locations = LetterPerPosition_Total %>% filter(Chr == '-') %>% group_by(Position) %>% 
                                                          mutate(Total = sum(Total)) %>% ungroup() %>% 
                                                          select(Position, Total) %>% distinct()
        
          #-Insert locations processing
        
        Insert_data = as.data.frame(indel(datos$aln)@insertion)
        DFGraph_Data <- as.data.frame(Graph_Data)
        Insert_data <- Insert_data %>% mutate(Count = DFGraph_Data['Abundance'][group,]) %>% select(start, end, Count)
        Insert_per_loci <- as.data.frame(c(1:max(width(datos$Sequences))))
        colnames(Insert_per_loci) <- c('Position')
        Insert_per_loci <- Unravel_Positions(Insert_per_loci, Insert_data)
        
        #-----------------Plotting--------------------------
        plot.new()
        
        #-Base Frequency
        
        colorscheme = c('red', 'black', 'blue', 'green', 'purple')
        colorscheme = setNames(colorscheme, c('T', 'G', 'C', 'A', '-'))
        
        MF <- plot_ly(LetterPerPosition_Norm, x = ~Position, y = ~Total, type = 'bar', name = ~Chr, color= ~Chr , colors = colorscheme) %>%
          layout(yaxis = list(title = 'Count'), barmode = 'stack')
        output$Mut_Freq <- renderPlotly({print(MF)}) 
        
        charts$tmpFileMF <- tempfile(fileext = '.pdf')
        export(MF, charts$tmpFileMF)
        
        #-Deletions per loci
        
        DL <- plot_ly(Deletions_locations, x = ~Position, y = ~Total, type = 'bar', name = 'Deletions') %>%
          layout(title = 'Deletions per position' , yaxis = list(title = 'Count'))
        output$Del_loc <- renderPlotly({print(DL)})
        
        charts$tmpFileDL <- tempfile(fileext = '.pdf')
        export(DL, charts$tmpFileDL)
        
        #-Deletion Sizes
        
        Del_Data = as.data.frame(InDel_Extraction(datos$aln, VCAbundance)[1])
        
        DS = plot_ly(Del_Data, y = ~TotalD, x = ~Deletions, type = 'bar' ) %>%
          layout(title = 'Deletions per sizes' , yaxis = list(title = 'Count'), xaxis = list(title = 'Size'))
        output$Del_sizes <- renderPlotly({print(DS)})
        
        charts$tmpFileDS <- tempfile(fileext = '.pdf')
        export(DS, charts$tmpFileDS)
        
        #-Insertion Sizes
        
        In_Data = as.data.frame(InDel_Extraction(datos$aln, VCAbundance)[2])
        
        IS = plot_ly(In_Data, y = ~TotalI, x = ~Insertions, type = 'bar' ) %>%
          layout(title = 'Insertions per sizes',yaxis = list(title = 'Count'), xaxis = list(title = 'Size')) 
        output$In_sizes <- renderPlotly({print(IS)})
        
        charts$tmpFileIS <- tempfile(fileext = '.pdf')
        export(IS, charts$tmpFileIS)
        
        #-Insertion per Loci
        
        IL = plot_ly(Insert_per_loci, y = ~Total_Count, x = ~Position, type = 'bar' ) %>%
          layout(title = 'Insertions per position', yaxis = list(title = 'Count'), xaxis = list(range = c(1, 250)))
        output$In_loc <- renderPlotly({print(IL)})
        
        charts$tmpFileIL <- tempfile(fileext = '.pdf')
        export(IL, charts$tmpFileIL)
        
        #- Dumped count
        
        output$Dump = renderPrint(paste('Filtered Clusters because of size restrictions: ', as.character(count(Dump)),sep = '')) #
    })
  })
  
  rows_selected <- reactive({ if( length(input$tabla2_rows_selected) == 0 ){  
    
        return('default')
    
    }else{ 
    
        return(input$tabla2_rows_selected )
    }})
  
  input_filename <- reactive({paste(str_replace(#Removes posible double underscore "_"
                              str_replace(#Removes Read number ID
                                substring(text = input$R1$name, 
                                          first = 1, 
                                          last = nchar(input$R1$name)-6)
                                ,'R1',''), '__',''), '_Results', sep='')})
  
  datos$alignmentdir <- 'empty' #Temporal asignment
  
  # datos$alignments <- eventReactive(input$tabla_rows_selected,{
  #   
  #   x <- list()
  #   cluster <- paste('Cluster', input$tabla_rows_selected, sep = '')
  #   
  #   if (!is.null(x$cluster)){
  #     x[[cluster]] <- datos$text
  #     
  #     return(x)
  #     }
  #   return(x)
  #   })
  

  output$downloadReport <- downloadHandler(
    filename = function() {
      paste(input_filename(),'pdf', sep = '.')
    },
    content = function(file) {
      src <- normalizePath('reports.Rmd')
      out <- rmarkdown::render( 'reports.Rmd', params = list(selection = rows_selected(),
                                                             table = datos$Tabla,
                                                             chartPie = charts$tmpFilePie,
                                                             chartDS = charts$tmpFileDS,
                                                             chartDL = charts$tmpFileDL,
                                                             chartIS = charts$tmpFileIS,
                                                             chartIL = charts$tmpFileIL,
                                                             chartMF = charts$tmpFileMF,
                                                             alignments = datos$alignmentdir) )
      file.rename(out, file) }
    )
  
  callModule(downloadCSV, 'downloadFile2', datos$Tabla, input_filename())
}

# Run the application
shinyApp(ui = ui, server = server)

