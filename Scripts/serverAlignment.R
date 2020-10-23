# server-Alignment calculation and printing
# Isolated script for the loading of the Alignment operations

#tabla_rows_selected() allow alignment selection regardless of using table 
#calculated with Reference or Target
tabla_rows_selected <- eventReactive(c(input$tablaR_rows_selected, 
                                       input$tablaT_rows_selected),{
                                        
                                        if(input$Alignment == 'Reference'){ 
                                         return( input$tablaR_rows_selected)
                                        }else if(input$Alignment == 'Target'){
                                         return( input$tablaT_rows_selected)
                                        }})
#----Single-Row Alignments----
datos$prev_selected_row <- 'none'
datos$prev_alnTo <- 'none'
datos$prev_alnType <- 'none'

alignment <- reactive({datos$text})

observeEvent( c( input$tablaR_rows_selected, input$tablaT_rows_selected, 
                 input$alnType, input$alnTo), {
                  
  req( input$tablaR_rows_selected != datos$prev_selected_row || input$alnType != datos$prev_alnType || input$alnTo != datos$prev_alnTo)
  
  withProgress( message = "Performing analisis", {
   
   req( datos$Tabla )
   req( tabla_rows_selected( ) )
   
   
   datos$prev_selected_row <- tabla_rows_selected()
   
   # tabla_rows_selected() provides the row that has been clicked, in the case of the unfiltered Tabla() there si no issue,
   # however if the user filters the data Tabla(), the row clicked does not correspond anymore with the Cluster#
   # Thus we need to find the ID in datos$Tabla_Original and then extract the sequence from datos$namedStringSet
   
   id_position <- names(datos$namedStringSet) %in% paste(Tabla()[tabla_rows_selected(),1],'_',rownames(Tabla())[tabla_rows_selected()],sep = '')
   Seq <- datos$namedStringSet[id_position]
   
   if(input$alnTo == "Reference")
   {
    datos$aln2 = pairwiseAlignment( datos$Ref, Seq, substitutionMatrix = nucleotideSubstitutionMatrix(match = 5, mismatch = -4),
                                    type = input$alnType, gapOpening=10, 
                                    gapExtension=0.5 )
   }else{
    datos$aln2 = pairwiseAlignment( datos$Target, Seq, substitutionMatrix = nucleotideSubstitutionMatrix(match = 5, mismatch = -4),
                                    type = input$alnType, gapOpening=10, 
                                    gapExtension=0.5)
   }
   
   diff = character(length( datos$aln2@pattern ))
   comp = character(length( datos$aln2@pattern ))
   ref = unlist(strsplit(as.character( datos$aln2@pattern ), 
                         split = "" ))
   query = unlist(strsplit(as.character( datos$aln2@subject ), 
                           split = "" ))
   
   for (i in 1:length( ref ))
   {
    if (ref[i] != query[i])
    {
     diff[i] = query[i]
     comp[i] = " "
    } else {
     diff[i] = "."
     comp[i] = "|"
    }
   }
   
   pos <- Position_vector(as.character(datos$aln2@pattern), 
                          ref) #Generates positioning vector
   
   if (!is.null(tabla_rows_selected()))
   {
    datos$text = as.data.frame(c(

     paste(pos, sep = '', collapse = ''), 
     paste(ref, sep = "", collapse = ""),
     paste(comp, sep = "", collapse = ""),
     paste(query, sep = "", collapse = ""),
     paste(diff, sep = "", collapse = "")

    ))
    
    if(input$alnTo == "Reference")
    {
     rownames(datos$text) = c("Position", "Reference", "Comparisson",
                              datos$Tabla$ID[tabla_rows_selected()], 
                              "Difference")
    }else{
     rownames(datos$text) = c("Position", "Target", "Comparisson",
                              datos$Tabla$ID[tabla_rows_selected()], 
                              "Difference")
    }
    
    colnames(datos$text) = c("Alignment")
    output$alignment = renderPrint(alignment())
    output$score = renderText(
     paste(
      "Alignment Score:",
      as.numeric(datos$aln2@score),
      "Alignment Length",
      length(ref),
      sep = "  ",
      collapse = ""
     )
    )
   }
   
   datos$clustering4StringSet <- datos$clustering %>% unnest() %>% filter(ClusterN == tabla_rows_selected()) %>% select(-ClusterN)
   datos$clustersStringSet <- DNAStringSet(datos$clustering4StringSet$SEQ)
   names(datos$clustersStringSet) <- datos$clustering4StringSet$ID
   
   output$Display_Unaligned <- renderUI({
    checkboxInput("display_unaligned", 
                  p(strong("Display unaligned fasta")), 
                  value = FALSE)
   })
   
   output$Blast_Search <- renderUI({column(3, wellPanel(
    helpText(a(div(img(src="BLASTn.png")), 
               href = paste('https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch',
                            '&JOB_TITLE=Cluster%20',
                            tabla_rows_selected(),
                            '&QUERY=>',
                            datos$Tabla$ID[tabla_rows_selected()],
                            '%0A',
                            Seq, sep=''), 
               target="_blank"))))
   })
   
   
   output$downloadFASTA_file <- renderUI({
    br()
    downloadButtonModule('downloadFASTA1',
                         'Download FASTA')
   })
   
   # output$downloadFASTAS_clstr <- renderUI({
   #  br()
   #  downloadButtonModule('downloadFASTAS_clstr',
   #                       'Download all Cluster FASTAs')
   # })
   
   br()
   output$downloadPAIR_file <- renderUI({
    br()
    downloadButtonModule('downloadPAIR1',
                         'Download Pairwise alignment')
   })
   
   output$Unaligned_Fasta <- renderUI({
    br()
    wellPanel(
     h5('Unaligned Sequence'),
     verbatimTextOutput('unalignedID'),
     verbatimTextOutput('unalignedFASTA'))
   })
   # output$downloadMSA <- renderUI({
   #  br()
   #  downloadButton("downloadMSA_clstr", 
   #                 "Download Cluster Multiple Sequence Alignment")
   #  
   # })
   
   # if(length(datos$clustersStringSet) > 100){
   # output$warningClustrUI <- renderUI(conditionalPanel(condition = "warningCLSTR() != FALSE",
   #  span(textOutput('warningCLSTR'), style="color:red"),
   #                                    br()))
   # }
   
   # output$warningCLSTR <- renderText(warningCLSTR())
   
   output$unalignedID <- renderText(as.character(
    datos$namedSequences[tabla_rows_selected(),]$ID))
   
   output$unalignedFASTA <- renderText(as.character(
    datos$namedSequences[tabla_rows_selected(),]$Sequence))
  })
  })