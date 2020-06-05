#
# 
# TO DO
# *Implementar una opcion para filtrar los primers. Estoy hay que hacerlo en la parte del prinseq. Los primers se meteran de forma manual
# *Indicar cual es la secuancia mas parecida al target.
# *Crear graficas de: * crispresso
#       - % de mutacion por base
#       - % de indel por base
#       - % cambio por base
#   * Crear un report en PDF (con markdown) con las graficas, la tabla y los primeros 10 alineamientos
#   * Trimm primer sequences (optional) - por numero de bases o por secuencia
#     -Error with certain datasets, fix ASAP
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#---Libraries----

library(shiny)
library(plotly)
library(tidyverse)

if (require('rmarkdown') == FALSE){
  
  install.packages('rmarkdown')
  library(rmarkdown)
  
}

#To install Biostrings:
# source("https://bioconductor.org/biocLite.R")
# biocLite("Biostrings")

library(Biostrings) 
library(DT)


#---UI---- 

ui <- fluidPage(# Application title
  titlePanel("CRISPRCas Analysis"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        h4("Input Files"),
        fileInput("R1", "Reads 1"),
        fileInput("R2", "Reads 2"),
        fileInput("Reference", "Reference", accept = c('.fasta', '.fastq','.txt')),
        fileInput("Target", "Target", accept = c('.fasta','.fastq','.txt')
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
          textInput('Adapter','Residual adapter trimming'),
          fileInput('Adapter', 'From file (Fasta format)', accept = '.fasta'),
          fluidRow(
            column(4,
                   actionButton('reset', 'Reset Input')
            ))
        ),
        numericInput('trimP1', 
                     '1st Primer trimmering',
                     min = 0,
                     value = 0,
                     max = 125),
        numericInput('trimP2', 
                     '2nd Primer trimmering',
                     min = 0,
                     value = 0,
                     max = 125)
      ),
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
          value = 0.9
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
    )),
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
                       max = 250*1.10,
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
                   wellPanel(
                   downloadButton("Report", "Generate report")
                   )))
                  )

  ))) 

#---SERVER----
# Define server logic required to draw a histogram
server <- function(input, output) {

  #---Reactive values----  
  options(shiny.maxRequestSize = 100 * 1024 ^ 2) #previous Mb value = 30
  datos <- reactiveValues()
  zoom <- reactiveValues(x = NULL, y = NULL) # Reactive values for zooming in on graphs
  values <- reactiveValues(upload_state = NULL)
  
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
  
  
  
  #------------EVENTS---------------------
  
  observeEvent(input$file1, {
      values$upload_state <- 'uploaded'
    })
  observeEvent(input$reset, {
      values$upload_state <- 'reset'
    })
    file_input <- reactive({
      if (is.null(values$upload_state)) {
        return(NULL)
      } else if (values$upload_state == 'uploaded') {
        return(input$file1)
      } else if (values$upload_state == 'reset') {
        return(NULL)
      }
    })
  
  observeEvent(input$Accept, {
    
    req(input$R1, input$R2, input$Reference,  input$Target) #Prevents App from crashing #Will work without target? no
    
    withProgress(message = "Performing analisis", {
      
      command = paste(
        CompletePATH,
        "/pipeline_v3.pl --r1 ",   ### Automatizar el PATH del script -- Doned
        input$R1$datapath,
        " --r2 ",
        input$R2$datapath,
        " --min_len ",
        input$MinLength,
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
      
      #Lets try reactive functions: not working
      
      # fasta <- reactive({
      #   validate(
      #    need(read_tsv("cluster.tsv", col_names = FALSE) != "", "All reads have been filtered, please loosen up restriction"))
      #   read_tsv("cluster.tsv")
      # })
      # 
      # datos$fasta = get(fasta())
      
      datos$fasta = read_tsv("cluster.tsv", col_names = FALSE)
      datos$cluster = read_tsv("cluster.bak.tsv", col_names = FALSE)
      
      
      if (isEmpty(colnames(datos$fasta))){#-If dataset is completely filtered due to any restriction app crashes
        stopApp(returnValue = 'All Reads have been filtered')
        }else{
        colnames(datos$fasta) = c("ID", "Sequence")
        colnames(datos$cluster) = c("ClusterN", "Length", "ID", "Identity")
        } #WE NEED TO BE ABLE TO SUPPLY ERROR INFORMATION ASAP
      
      #colnames(datos$fasta) = c("ID", "Sequence")
      #colnames(datos$cluster) = c("ClusterN", "Length", "ID", "Identity")
      
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
      
      system("rm join* good* bad* cluster*"); #To prevent constant cashing on the same data
      
      datos$Tabla = inner_join(datos$Tabla,tmp) %>% mutate(score = round(score,1), Freq = signif(Freq,2)) %>% arrange(desc(Abundance))

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
      
      for (i in 1:length(ref))
      {
        if (ref[i] != query[i])
        {
          comp[i] = query[i]
        } else {
          comp[i] = "."
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
      output$Blast_Search <- renderUI({a("Blast Search", 
                                         href = paste('https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch',
                                                      '&JOB_TITLE=Cluster%20',
                                                      input$tabla_rows_selected,
                                                      '&QUERY=>',
                                                      datos$Tabla$ID[input$tabla_rows_selected],
                                                      '%0A',
                                                      Seq, sep=''), 
                                         target="_blank")})
    })

  })
  observeEvent(c(input$Report),{
               
      output$report <- downloadHandler(
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

