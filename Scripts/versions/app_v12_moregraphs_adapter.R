#
# 
# TO DO
# *Indicar cual es la secuancia mas parecida al target.
# *Crear graficas de: * crispresso
#       - % de mutacion por base
#       - % de indel por base
#       - % cambio por base
#   * Crear un report en PDF (con markdown) con las graficas, la tabla y los primeros 10 alineamientos
# 
# *Check how the cluster is determining the consensus sequence
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
source(paste(CompletePATH, '/../Modules/downloadCSVModule.R', sep = ''))


#---UI Side Panel---- 

ui <- fluidPage(# Application title
  
  titlePanel( title=(div("CRISPRCas Analysis", img(src="dna_free.png")))),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        h4("Input Files"),
        fileInput("R1", "Reads 1", accept = c('.fastq', '.fastq.gz')),
        fileInput("R2", "Reads 2", accept = c('.fastq', '.fastq.gz')),
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
                                "5' sequence",
                                min = 0,
                                value = 0,
                                max = 125),
                   numericInput('trimA2', 
                                "3' sequence",
                                min = 0,
                                value = 0,
                                max = 125)),
          tabPanel(p(strong('by sequence')),
            tabsetPanel(id = 'input_type_A',
              tabPanel('Text input', value = 'Direct_input_A',
                    br(),
                    textInput("R1_Adapter","R1 sequence"),
                    textInput("R2_Adapter","R2 sequence"),
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
            verbatimTextOutput(outputId = 'checkR2A'))
          
        )),
        wellPanel(
        navbarPage(title = p(strong('Primer trimming:')),
          tabPanel(p(strong('None'))),
          tabPanel(p(strong("by length")),
                   numericInput('trimP1', 
                                         "5' sequence",
                                        min = 0,
                                        value = 0,
                                        max = 125),
                   numericInput('trimP2', 
                                        "3' sequence",
                                        min = 0,
                                        value = 0,
                                        max = 125)),
          tabPanel(p(strong(("by sequence"))),
                   tabsetPanel(id = 'input_type_P',
                     tabPanel("Text input",
                              value = 'Direct_input_P',
                        hr(),
                        textInput('F_Primer1',"5' forward end sequence"),
                        textInput('R_Primer2',"5' reverse end sequence"),
                        br(),
                        fluidRow(#column(4, actionButton('resetP', 'Reset Input')), #No need for this anymore
                                 column(4, actionButton('uploadPrimers', 'Upload Primers'))),
                        br(),
                        br()
                        ),
                     tabPanel('File input', value = 'File_input_P',
                        hr(),
                        fileInput('PrimerFile', 'From file (Fasta format)', accept = '.fasta'),
                        p(em('The file should contain four lines, headers on line 1 and 3 starting with an arrow ">", 
                              first sequence being the 5" primer for the Forward strand and second sequence the 5" primer for the 
                              Reverse. The rest of the reads are
                              obtained automatically.')),
                        br(),
                        br())),
                        verbatimTextOutput(outputId = 'checkFP1'),
                        verbatimTextOutput(outputId = 'checkFP2'),
                        verbatimTextOutput(outputId = 'checkRP1'),
                        verbatimTextOutput(outputId = 'checkRP2')
            
          
        )
      ))),
  #---Well stablished parameters----
      wellPanel(
        h4("FastQ-Join parameters"),
        
        numericInput(
          "Np",
          "Maximun difference (%)",
          min = 1,
          max = 100,
          value = 0
        ),
        numericInput(
          "Nm",
          "Minimun overlap (nucleotides)",
          min = 4,
          max = 999,
          value = 10
        )
      ),
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
      
      actionButton("Accept", "Compute")
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
                    verbatimTextOutput("alignment"),
                    verbatimTextOutput("score"),
                    uiOutput('Blast_Search')
                   
                   ),
                   
          tabPanel("Graphics", value = 'Graphics',  # Show a plot of the generated distribution
                    br(),
                    wellPanel(
                      actionButton(inputId =  "GG", "Generate Plots")
                    ),
                   
                    wellPanel(
                      h4('Mutation Frequencies'),
                      plotlyOutput(outputId = 'MutFreq')
                    ),
                   
                    wellPanel(
                      h4('Deletions per position'),
                      plotlyOutput(outputId = 'Del_loc')
                    ),
                  
                    wellPanel(
                      h4('Deletions sizes'),
                      plotlyOutput(outputId = 'Del_sizes')
                    ),
                   
                   wellPanel(
                     h4('Insertion sizes'),
                     plotlyOutput(outputId = 'In_sizes')
                   ),
                   
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
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  datos <- reactiveValues()
  file_upload <- reactiveValues(upload_Adapter_state = NULL)
  file_upload <- reactiveValues(upload_Primer_state = NULL)
  R1_file <- reactive({
    if (str_sub(input$R1$datapath, -3) == '.gz'){
      Unzip_file(CompletePATH, input$R1$datapath)}
    else{
      input$R1$datapath}
    })
  R2_file <- reactive({
    if(str_sub(input$R1$datapath, -3) == '.gz'){
      Unzip_file(CompletePATH, input$R1$datapath)}
    else{
      input$R1$datapath}
  })

  #----Filtering Parameters declarations----
  Target_FileInput <- callModule(inputReset_server, 'reset_Target')
  #-Prevents the appereance of the Target alignment option if no Target has been uploaded
  output$is.nullTarget <- renderUI({#Generates reactive UI depending on the presence of a Target or not. / not working
    selectInput(
      "alnTo",
      "Alignment to:",
      if (!is.null(Target_FileInput())){c("Reference","Target")}else{c("Reference")})
  })
  
  #This following set of observeEvents is to allow user to change his mind on uploaded Adapter files, 
  #these files currently are not going anywhere though.
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
  
  output$checkFP1 <- renderPrint(paste("Forward read 5': ",PrimerR1_5(), sep = ''))
  output$checkFP2 <- renderPrint(paste("Forward read 3': ",PrimerR1_3(), sep = ''))
  output$checkRP1 <- renderPrint(paste("Reverse read 5': ",PrimerR2_5(), sep = ''))
  output$checkRP2 <- renderPrint(paste("Reverse read 3': ",PrimerR2_3(), sep = ''))

  #ToDo: Introduce a check to make sure length of primers is not over 30?
  #----Processing----
  
  observeEvent(input$Accept, {
    
    req(input$R1, input$R2, input$Reference)
    
    
    withProgress(message = "Performing analisis", {
      
      command = paste(
        CompletePATH,
        "/pipeline_v5.pl --r1 ",   ### Automatizar el PATH del script -- Doned
        R1_file(),
        " --r2 ",
        R2_file(),
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
      
      datos$Ref = readDNAStringSet(input$Reference$datapath)
      
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
      
      system("rm R1_* R2_* join* good* bad* cluster*"); #To prevent constant cashing on the same data
      
      datos$Tabla = inner_join(datos$Tabla_raw,tmp) %>% mutate(score = round(score,1), Freq = signif(Freq,2)) %>% 
        arrange(desc(Abundance))

        output$tabla = renderDT(datos$Tabla, selection = "single")
        output$tabla2 = renderDT(datos$Tabla) #Table with multiple selection for PDF formatting
        output$Print = renderPrint(datos$Tabla[sort.default(input$tabla2_rows_selected),])
        
      
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
        text = as.data.frame(c(
          paste(pos, sep = '', collapse = ''), #Not working well enough, fails in presence of inserts
          paste(ref, sep = "", collapse = ""),
          paste(query, sep = "", collapse = ""),
          paste(comp, sep = "", collapse = "")
        ))
        
        if(input$alnTo == "Reference")
        {
          rownames(text) = c("Position", "Reference", datos$Tabla$ID[input$tabla_rows_selected], "Comparisson")
        }else{
          rownames(text) = c("Position", "Target", datos$Tabla$ID[input$tabla_rows_selected], "Comparisson")
        }
        
        colnames(text) = c("Alignment")
        output$alignment = renderPrint(text)
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
    
          #Necesary for Indel data extraction
        VCAbundance = select(datos$Tabla, Abundance)
        
        Insert_data = as.data.frame(indel(datos$aln)@insertion)
        unordered_ID = as.data.frame(names(datos$Sequences))
        colnames(unordered_ID) <- 'names'
        Identifiers = unordered_ID %>% separate(names, c('names', 'kk'), sep = ' ') %>% select(-kk)
          #Now contrast with datos$Tabla to obtain Abundance
        
        
        
          #Extracting Sequences, identificators, sizes
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
                                      
        
        LetterPerPosition_Norm = LetterPerPosition_Total %>% ungroup() %>% mutate_at(vars(Total),funs(./sum(Total))) %>% 
                                      ungroup(Position) %>% group_by(Chr, Position) %>%
                                      mutate_at(vars(Total), funs(sum(Total))) %>% 
                                      group_by(Total, Position) %>% distinct(Chr)
        
        #-Isolating deletions
        Deletions_locations = LetterPerPosition_Total %>% filter(Chr == '-')
        
        #-----------------Plotting--------------------------
        plot.new()
        
        #-Base Frequency
        
        MF <- plot_ly(LetterPerPosition_Norm, x = ~Position, y = ~Total, type = 'bar', name = ~Chr) %>%
          layout(yaxis = list(title = 'Count'), barmode = 'stack')
        output$MutFreq <- renderPlotly({print(MF)}) 
        
        #-Deletions per loci
        
        DL <- plot_ly(Deletions_locations, x = ~Position, y = ~Total, type = 'bar', name = 'Deletions') %>%
          layout(yaxis = list(title = 'Count'))
        output$Del_loc <- renderPlotly({print(DL)})
        
        #-Deletion Sizes
        
        Del_Data = as.data.frame(InDel_Extraction(datos$aln, VCAbundance)[1])
        
        DS = plot_ly(Del_Data, y = ~TotalD, x = ~Deletions, type = 'bar' )
        output$Del_sizes <- renderPlotly({print(DS)})
        
        #-Insertion Sizes
        
        In_Data = as.data.frame(InDel_Extraction(datos$aln, VCAbundance)[2])
        
        IS = plot_ly(In_Data, y = ~TotalI, x = ~Insertions, type = 'bar' )
        output$In_sizes <- renderPlotly({print(IS)})
        
        #-Insertion per Loci
        
        output$Dump = renderPrint(paste('Filtered Clusters because of size restrictions: ', as.character(count(Dump)),sep = '')) #
  })
  

  output$downloadReport <- downloadHandler(
    filename = function() {
      paste('my-report','pdf', sep = '.')
    },
    content = function(file) {
      src <- normalizePath('reports.Rmd')
      out <- rmarkdown::render( 'reports.Rmd', pdf_document(), params = list(selection = input$tabla2_rows_selected, table = datos$Tabla) )
      file.rename(out, file) }
    )
  filename <- reactive({paste(substring(text = input$R1$name, first = 1, last = nchar(input$R1$name)-6), '_Results', sep='')})
  callModule(downloadCSV, 'downloadFile2', datos$Tabla, filename())
}

# Run the application
shinyApp(ui = ui, server = server)

