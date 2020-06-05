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

library(shiny)
library(plotly)
library(tidyverse)

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



#---UI Side Panel---- 

ui <- fluidPage(# Application title
  
  titlePanel( title=(div("CRISPRCas Analysis", img(src="dna_free.png")))),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        h4("Input Files"),
        fileInput("R1", "Reads 1", accept = '.fastq'),
        fileInput("R2", "Reads 2", accept = '.fastq'),
        fileInput("Reference", "Reference", accept = c('.fasta', '.fastq')),
        fileInput("Target", "Target", accept = c('.fasta','.fastq'))
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
            tabsetPanel(id = 'Adapter by sequence',
              tabPanel('Text input',
                   br(),
                   textInput("Adapter 5'","5' sequence"),
                   textInput("Adapter 3'","3' sequence")),
              tabPanel('File input',
                   fileInput("Adapter", 'From file (Fasta format)', accept = '.fasta'),
                   fluidRow(column(4, actionButton('resetA', 'Reset Input')))
                   )),
            
            verbatimTextOutput(outputId = 'checkFA1'),
            verbatimTextOutput(outputId = 'checkFA2'),
            verbatimTextOutput(outputId = 'checkRA1'),
            verbatimTextOutput(outputId = 'checkRA2'))
          
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
                   tabsetPanel(id = 'input_type',
                     tabPanel("Text input",
                              value = 'Direct_input',
                        hr(),
                        textInput('F_Primer1',"5' forward end sequence"),
                        textInput('R_Primer2',"5' reverse end sequence"),
                        br(),
                        fluidRow(#column(4, actionButton('resetP', 'Reset Input')), #No need for this anymore
                                 column(4, actionButton('uploadPrimers', 'Upload Primers'))),
                        br(),
                        br()
                        ),
                     tabPanel('File input', value = 'File_input',
                        hr(),
                        fileInput('PrimerFile', 'From file (Fasta format)', accept = '.fasta'),
                        p(em('The file should contain four lines, headers on line 1 and 3 starting with an arrow ">", 
                              first sequence being the 5" primer for the Forward strand and second sequence the 5" primer for the 
                              Reverse. The rest of the reads are
                              obtained automatically.')),
                        br(),
                        #fluidRow(column(4, actionButton('resetP', 'Reset Input'))), #No need for this anymore
                        #br(),
                        #br(),
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
          value = 8
        ),
        numericInput(
          "Nm",
          "Minimun overlap (nucleotides)",
          min = 4,
          max = 200,
          value = 6
        )
      ),
      wellPanel(
        h4("Clustering Parameter"),
        sliderInput(
          "cov",
          "Coverage",
          min = 0 ,
          max = 1,
          value = 0.99
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
    
        tabsetPanel(
          tabPanel("Clusterization & Alignment",
                   
                   h3("Table"),
                   DTOutput("tabla"),
                   h3("Alignment"),
                   selectInput(
                     "alnType",
                     "Alignment method",
                     c("global", "local", "overlap", "global-local", "local-global")
                   ),
                   selectInput(
                     "alnTo",
                     "Alignment to:",
                     c("Reference","Target")
                     
                   ),
                    verbatimTextOutput("alignment"),
                    verbatimTextOutput("score"),
                    uiOutput('Blast_Search')
                   
                   ),
                   
          tabPanel("Graphics",  # Show a plot of the generated distribution
                    br(),
                    wellPanel(
                     actionButton(inputId =  "GG", "Generate Mutation Frequency Plot")
                     ),
                   
                    wellPanel(
                    plotlyOutput(outputId = 'MutFreq')
                     ),
                    column(width = 4,
                      h4("Zoomin parameters"),
                         
                      numericInput(
                       "X1",
                       "Starting position",
                       min = 1,
                       max = 250*1.10, #Do we still need this?
                       value = 1 ),
                      numericInput(
                       'X2',
                       "End position",
                       min = 1,
                       max = 250*1.10,
                       value = 250),
                      textOutput(p(em(outputId = 'Dump'))),
                      textOutput(outputId = 'Dump')
                     )),
          tabPanel('Download',
                   column(width = 5,
                    br(),
                    wellPanel(
                    downloadButton("Report", "Generate report")
                   )))
                  )
  )))

#---SERVER---- 
# Define server logic required to draw a histogram
server <- function(input, output) {

  #---Reactive values----  
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  datos <- reactiveValues()
  zoom <- reactiveValues(x = NULL, y = NULL) # Reactive values for zooming in on graphs
  file_upload <- reactiveValues(upload_Adapter_state = NULL)
  file_upload <- reactiveValues(upload_Primer_state = NULL)
  

  #------Setting up work environment-----
  
  appPATH = getwd()
  
  if (substring(appPATH, nchar(appPATH)-8, nchar(appPATH)) != 'CrisPRAL'){
    
    system('find /home/ -name "CrisPRAL" > PATH.txt')
    
    appPATH = read_file('PATH.txt')
    
    system('rm PATH.txt')
    
    # ContainingFolder = substring(appPATH, 1, nchar(appPATH)-9) #-----14/09 error loading the correct pipeline, trying to switch to CrisPRAL directory
    # DataStorage = "analisisNGSCRISPRCas/"
    # CompletePATH = paste(ContainingFolder, DataStorage, sep = "")
    
    CompletePATH = substring(appPATH, 1, nchar(appPATH)-1) # nchar - 1 to remove "/n" due to .txt parsing
    
    setwd(CompletePATH)
    
  }else{CompletePATH = appPATH}
  
  
  
  
  #----Filtering Parameters declarations----
  
  #This following set of observeEvents is to allow user to change his mind on uploaded Adapter files, 
  #these files currently are not going anywhere though.
  observeEvent(input$Adapter, {
    file_upload$upload_Adapter_state <- 'uploaded'
    }) 
  observeEvent(input$resetA, {
    file_upload$upload_Adapter_state <- 'reset'
    })
  Adapter_FileInput <- reactive({
      if (is.null(file_upload$upload_Adapter_state)) {
        return(NULL)
      } else if (file_upload$upload_Adapter_state == 'uploaded') {
        return(input$Adapter)
      } else if (file_upload$upload_Adapter_state == 'reset') {
        return(NULL)
      }
    })
  #Same than the Adapter strings but with the Primers.
  observeEvent(c(input$PrimerFile, input$uploadPrimers), {
    file_upload$upload_Primer_state <- 'uploaded'
  }) 
  observeEvent(c(input$resetP, input$input_type), {
    file_upload$upload_Primer_state <- 'reset'
  })
  
  Primer_Directinput <- eventReactive(input$uploadPrimers,{ c(input$F_Primer1,input$R_Primer2) })
  Primer_FileInput <- reactive({
    if (is.null(file_upload$upload_Primer_state)) {
      return(NULL)
    } else if (file_upload$upload_Primer_state == 'uploaded') {
      
      if (input$input_type == 'File_input'){return(input$PrimerFile)}
      else{return(Primer_Directinput())}
        }
     else if (file_upload$upload_Primer_state == 'reset') {
      return(NULL)
    }
  })
  
  
  #Primer upload through file:
  #Primer sequences are loaded as reactive functions to prevent residual data lingering.
  Primer5_sequences <- reactive({
    #validate(need(!is.null(Primer_FileInput()),'Valid primers have not been uploaded'))
    if (is.null(Primer_FileInput())){
      return("Empty")
    }else if (!is.null(Primer_FileInput())){
      if(input$input_type == 'File_input'){return(readDNAStringSet(Primer_FileInput()$datapath))
      }else{
          return(DNAStringSet(Primer_FileInput()))}
      }})
  Primer3_sequences <- reactive({
    #validate(need(!is.null(Primer_FileInput()),'Valid primers have not been uploaded'))
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
    
    req(input$R1, input$R2, input$Reference,  input$Target) #Prevents App from crashing #Will work without target? no
    
    withProgress(message = "Performing analisis", {
      
      command = paste(
        CompletePATH,
        "/pipeline_v4.pl --r1 ",   ### Automatizar el PATH del script -- Doned
        input$R1$datapath,
        " --r2 ",
        input$R2$datapath,
        " --min_len ",
        input$MinLength,
        " --F_Primer1 ",
        as.character(PrimerR1_5()),
        " --R_Primer2 ",
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
      
      system(command)
      
      datos$fasta = read_tsv("cluster.tsv", col_names = FALSE)
      datos$cluster = read_tsv("cluster.bak.tsv", col_names = FALSE)
      
      if (isEmpty(colnames(datos$fasta))){#-If dataset is completely filtered due to any restriction app crashes
        stopApp(returnValue = 'Pipe was broken')
        }else{
        colnames(datos$fasta) = c("ID", "Sequence")
        colnames(datos$cluster) = c("ClusterN", "Length", "ID", "Identity")
        } #WE NEED TO BE ABLE TO SUPPLY ERROR INFORMATION ASAP
      
      datos$Tabla = datos$cluster %>% 
        group_by(ClusterN) %>% 
        mutate(Abundance = n()) %>% 
        ungroup() %>% 
        mutate(Freq = 100 *Abundance / n()) %>%
        filter(Identity == "*") %>% 
        select(ID, Abundance, Freq) %>% 
        arrange(desc(Abundance))
      
      datos$Ref = readDNAStringSet(input$Reference$datapath)
      datos$Target = readDNAStringSet(input$Target$datapath)
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
      
      datos$Tabla = inner_join(datos$Tabla,tmp) %>% mutate(score = round(score,1), Freq = signif(Freq,2)) %>% 
        arrange(desc(Abundance))

      output$tabla = renderDT(datos$Tabla, selection = "single")
      
    })
  
  })
  observeEvent(c(input$GG), { #---Plotting Graphs

        req(datos$Tabla) # Prevents APP from crashing
    
        #---Zooming options
        zoom$x1 <- input$X1
        zoom$x2 <- input$X2
        
        
        #-----------Setting up the Data----------------
        
        VCSeq = list(as.character(datos$aln@pattern))
        VCIDs = select(datos$Tabla, ID)
        SizeList = list(nchar(datos$aln@pattern))
        
        VCIDs = bind_cols(VCIDs, SizeList, VCSeq) #Bind Identificator with sequence and size as to not lose them when we filter by size
        colnames(VCIDs) = c('ID', 'Size', 'Sequence') #Name the tables
        
        FilteredVCIDs = filter(VCIDs, Size < 250*(2-input$cov)) #Leave enough margin in order to correct for coverage lax
        Dump = filter(VCIDs, Size > 250*(2-input$cov))
        
        TableSize = max(FilteredVCIDs['Size'])
  
        #----------------Processing-------------------------
        
        IDsSize = (count(FilteredVCIDs))*TableSize
        
        LetterPerPosition = data.frame(ID=rep(FilteredVCIDs$ID,TableSize), 
                                       Position = sort(rep(1:TableSize,count(FilteredVCIDs))),
                                       Chr = substring(FilteredVCIDs$Sequence, 
                                                       sort(rep(1:TableSize,count(FilteredVCIDs))), 
                                                       sort(rep(1:TableSize,count(FilteredVCIDs)))))
        
        LetterPerPosition = LetterPerPosition %>% group_by(Position, Chr) %>% count(Chr)
        LetterPerPosition = filter(LetterPerPosition, Chr != "")
        LetterPerPosition = rename(LetterPerPosition, n = 'Count')
        LetterPerPosition2 = LetterPerPosition %>% group_by(Position) %>% mutate_at(vars(Count),funs(./sum(Count)))
        
        
        #-----------------Plotting--------------------------
        plot.new()
        
        p <- plot_ly(LetterPerPosition2, x = ~Position, y = ~Count, type = 'bar', name = ~Chr) %>%
          layout(yaxis = list(title = 'Count'), barmode = 'stack')
        
        output$MutFreq <- renderPlotly({
          print(p)
        })
        output$Dump = renderPrint(paste('Filtered reads because of size restrictions: ', as.character(count(Dump)),sep = '')) #
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
      pos = c()
      inserts = 0
      
      for (i in 1:length(ref))
      {
        if (ref[i] != query[i])
        {
          comp[i] = query[i]
          if (ref[i] == '-'){
            inserts = inserts+1
          }else{
            if (i == 1 || i%%20 == 0){
              pos[i] = i-inserts
            }else{
              pos[i] = ' '
            }
          }
        } else {
          comp[i] = "."
          if (i == 1 || i%%20 == 0){
            pos[i] = i-inserts
          }else{
            pos[i] = ' '
          }
        }
      }
      
      if (!is.null(input$tabla_rows_selected))
      {
        text = as.data.frame(c(
          paste(ref, sep = "", collapse = ""),
          paste(query, sep = "", collapse = ""),
          paste(comp, sep = "", collapse = "")
        ))
        
        if(input$alnTo == "Reference")
        {
          rownames(text) = c("Reference", datos$Tabla$ID[input$tabla_rows_selected], "Comparisson")
        }else{
          rownames(text) = c("Target", datos$Tabla$ID[input$tabla_rows_selected], "Comparisson")
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
  observeEvent(c(input$Report),{
               
      output$Report <- downloadHandler(
         # For PDF output, change this to "report.pdf"
         filename = "report.pdf",
         content = function(file){
           
           # Copying the report file to a temporary directory before processing it, in
           # case we don't have write permissions to the current working dir.
           
           tempReport <- file.path(tempdir(), "report.Rmd")
           file.copy("report.Rmd", tempReport, overwrite = TRUE)
           params = 'na'
           rmarkdown::render(tempReport, output_file = file,
                             params = params,
                             envir = new.env(parent = globalenv())
           ) 
           
         }
      )
  })
  
  }

# Run the application
shinyApp(ui = ui, server = server)

