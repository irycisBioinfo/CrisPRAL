#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(tidyverse)
library(Biostrings)
library(DT)


# Define UI for application that draws a histogram
ui <- fluidPage(# Application title
  titlePanel("CRISPRCas Analysis"),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        h4("Input Files"),
        fileInput("R1", "Reads 1"),
        fileInput("R2", "Reads 2"),
        fileInput("Reference", "Reference")
      ),
      wellPanel(
        h4("Filter Reads Options"),
        numericInput(
          "MinLength",
          "Minimum reads length",
          value = 200,
          min = 31,
          max = 400,
          step = 10
        )
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
      actionButton("Accept", "Accept")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      h3("Table"),
      DTOutput("tabla"),
      h3("Alignment"),
      selectInput(
        "alnType",
        "Alignment method",
        c("global", "local", "overlap", "global-local", "local-global")
      ),
      verbatimTextOutput("alignment"),
      verbatimTextOutput("score")
      
    )
  ))

# Define server logic required to draw a histogram
server <- function(input, output) {
  setwd("/storage/analisisNGSCRISPRCas/")
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  datos <- reactiveValues()
  
  observeEvent(input$Accept, {
    withProgress(message = "Performing analisis", {
      command = paste(
        "/storage/analisisNGSCRISPRCas/pipeline.pl --r1 ",
        input$R1$datapath,
        " --r2 ",
        input$R2$datapath,
        " --min_len ",
        input$MinLength,
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
      colnames(datos$fasta) = c("ID", "Sequence")
      colnames(datos$cluster) = c("ClusterN", "Length", "ID", "Identity")
      datos$Tabla = datos$cluster %>% 
        group_by(ClusterN) %>% 
        mutate(Abundance = n()) %>% 
        ungroup() %>% 
        mutate(Freq = 100 *Abundance / n()) %>%
        filter(Identity == "*") %>% 
        select(ID, Abundance, Freq) %>% 
        arrange(desc(Abundance))
      datos$Ref = readDNAStringSet(input$Reference$datapath)
      
      output$tabla = renderDT(datos$Tabla, selection = "single")
      
    })
  })
  observeEvent(c(input$tabla_rows_selected, input$alnType), {
    withProgress(message = "Performing analisis", {
      if (is.null(datos$Tabla) | is.null(input$tabla_rows_selected))
      {
        return(NULL)
      }
      S = datos$Tabla %>% filter(ID == datos$Tabla$ID[input$tabla_rows_selected]) %>% left_join(datos$fasta)
      Seq = DNAStringSet(S$Sequence)
      aln = pairwiseAlignment(datos$Ref, Seq, type = input$alnType)
      
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
        rownames(text) = c("Reference", datos$Tabla$ID[input$tabla_rows_selected], "Comparisson")
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
      
    })
  })
  
  
  
}

# Run the application
shinyApp(ui = ui, server = server)
