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
   req( datos$Tabla )
   req( tabla_rows_selected( ) )
   
   
   datos$prev_selected_row <- tabla_rows_selected()
   
   # tabla_rows_selected() provides the row that has been clicked, in the case of the unfiltered Tabla() there is no issue,
   # however if the user filters the data Tabla(), the row clicked does not correspond anymore with the Cluster#
   # Thus we need to find the ID in datos$Tabla_Original and then extract the sequence from datos$namedStringSet
   
   id_position <- names(datos$namedStringSet) %in% paste(Tabla()[tabla_rows_selected(),1],'_',rownames(Tabla())[tabla_rows_selected()],sep = '')
   Seq <- datos$namedStringSet[id_position]
   
   if(input$alnTo == "Reference")
   { Subject <- datos$Ref }else{ Subject <- datos$Target }
   
    datos$aln2 = pairwiseAlignment( Seq, Subject, substitutionMatrix = nucleotideSubstitutionMatrix(match = 5, mismatch = -4),
                                    type = input$alnType, gapOpening=input$gap_open, 
                                    gapExtension=input$gap_extend )

   alignment = Fix_gaps(datos$aln2, Seq, Subject)
   query <- alignment$query[[1]]
   ref <- alignment$ref[[1]] 
   
   
   diff = character(length( query ))
   comp = character(length( query ))
   
   for(i in c(1:length( ref ))){
     if (ref[i] != query[i])
     {
       diff[i] = query[i]
       comp[i] = " "
     } else {
       diff[i] = "."
       comp[i] = "|"
     }
   }
   pos <- Position_vector(paste(query, collapse = ''), paste(ref, collapse = ''), datos$aln2) #Generates positioning vector from strings.
   
   if (!is.null(tabla_rows_selected()))
   {
     datos$text = as.data.frame(c(
       
       paste(pos, sep = '', collapse = ''),
       paste(ref, sep = "", collapse = ""),
       paste(comp, sep = "", collapse = ""),
       paste(query, sep = "", collapse = ""),
       paste(diff, sep = "", collapse = "")
        
        # paste(alignment, sep = '')
       
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
   
   output$Display_Unaligned <- renderUI({
    checkboxInput("display_unaligned", 
                  p(strong("Display unaligned fasta")), 
                  value = FALSE)
   })
   
   # ############################################################################.
   ##### MISMATCHES TABLE ####
   ##############################################################################.
   
   datos$TablaMM <- tibble(SNP = str_split(Table_to_download()[tabla_rows_selected(),] %>% 
                                              pull(SNP), pattern = ';')[[1]]) %>% 
      separate(SNP, c('lead','Alternate'), sep = '>') %>% 
      mutate(Reference = str_extract(lead, pattern = "[ACTG]$")) %>% relocate(Reference, .before = Alternate) %>%
      mutate(Position = as.numeric(str_remove(lead, pattern = "[ATCG]$"))) %>% 
      select(-lead) %>% relocate(Position, .before = Reference) %>%
      mutate(Nomenclature = str_c(Position,Reference,'>',Alternate)) %>%
      mutate(Type = 'Variant')
   
   # ############################################################################.
   # INDELS TABLES ####
   ##############################################################################.
   
   # Find selected entry in datos$aln
   alignment_selected <- match(Table_to_download()[tabla_rows_selected(),]$ID, datos$Tabla_unsort$ID)
   original_pos <- Position_vector(paste(query, collapse = ''), paste(ref, collapse = ''), datos$aln[alignment_selected])
   
   #################### ROUGH NOTES & TRIALS #######################.
   
   del_pos <- original_pos
   ins_pos <- original_pos
   
   adjustment_d <- lag(cumsum(width(deletion(datos$aln[alignment_selected])[[1]])))
   adjustment_d[1] <- 0
   map_pos_dels_aln <- start(deletion(datos$aln[alignment_selected])[[1]])+adjustment_d
   
   alt_dels <- str_sub(datos$aln@pattern[alignment_selected],map_pos_dels_aln-1,map_pos_dels_aln-1)
   ref_dels <- str_sub(datos$aln@subject[alignment_selected],map_pos_dels_aln-1,map_pos_dels_aln+(width(deletion(datos$aln[alignment_selected])[[1]])-1))
   
   adjustment_i <- lag(cumsum(width(insertion(datos$aln[alignment_selected])[[1]])))
   adjustment_i[1] <- 0
   map_pos_ins_aln <- start(insertion(datos$aln[alignment_selected])[[1]])+adjustment_i
   
   alt_ins <- str_sub(datos$aln@subject[alignment_selected],map_pos_ins_aln-1,map_pos_ins_aln-1)
   ref_ins <- str_sub(datos$aln@pattern[alignment_selected],map_pos_ins_aln-1,map_pos_ins_aln+(width(insertion(datos$aln[alignment_selected])[[1]])-1))
   
   ##################################################.
   
   datos$TablaDels <- tibble(deletion_range = str_split(Table_to_download()[tabla_rows_selected(),] %>% 
                                                           select(deletion_range) %>% 
                                                           pull(), pattern = ';')[[1]]) %>% 
                        separate(deletion_range, c('Start','End'), sep = ':') %>% mutate(Size = c(as.numeric(End)-(as.numeric(Start)-1)))
   
   del_pos[!str_split(datos$aln@subject[alignment_selected], '')[[1]] %in% '-'][as.numeric(datos$TablaDels$Start) - 1] <- as.numeric(datos$TablaDels$Start) - 1
   del_pos_maped <- match(as.numeric(datos$TablaDels$Start) - 1,del_pos)
   
   datos$TablaDels <- datos$TablaDels %>% 
      mutate(Position = as.numeric(Start) -1) %>% 
      mutate(Reference = str_sub(unaligned(datos$aln@subject[alignment_selected]), start = Position, end = End)) %>% 
      mutate(Alternate = str_sub(unaligned(datos$aln@pattern[alignment_selected]), start = Position, end = Position)) %>% 
      select(-End, -Start) %>% mutate(Nomenclature = str_c(Position, Reference,'>',Alternate, sep = '')) %>%
      mutate(Type = 'Deletion')
   datos$TablaIns <- tibble(insert_range = str_split(Table_to_download()[tabla_rows_selected(),] %>% 
                                                           select(insert_range) %>% 
                                                           pull(), pattern = ';')[[1]]) %>% 
      separate(insert_range, c('Start','End'), sep = ':')
   
   datos$TablaIns <- datos$TablaIns %>% 
      mutate(Position = as.numeric(Start) -1) %>% 
      mutate(Reference = str_sub(unaligned(datos$aln@subject[alignment_selected]), start = Position, end = Position)) %>% 
      mutate(Alternate = str_sub(unaligned(datos$aln@pattern[alignment_selected]), start = Position, end = End)) %>% 
      select(-End,-Start) %>% mutate(Nomenclature = str_c(Position, Reference,'>',Alternate, sep = '')) %>%
      mutate(Type = 'Insert')
   
   ################# Fuse all data ##################.
   
   datos$variants <-  bind_rows(datos$TablaMM,datos$TablaDels,datos$TablaIns) %>% drop_na() %>% arrange(Position)
   
   ################ Define dinamic UI ###############.
   
   output$mutations <- renderUI({
      
      wellPanel(
         checkboxInput('show_indels',label = 'Display formatted mismatches and InDels', value = TRUE),
         conditionalPanel(condition = 'input.show_indels == true',
                          fluidRow(
                                    wellPanel(
                                       h3("Mutations"),
                                       DTOutput('formatted_variants'))
                          ) # fluidRow end bracket
                          
         ))} # renderUI end bracket
   )
   
   # ############################################################################.
   ##### BLAST SEARCH AND DOWNLOADS ####
   ##############################################################################.
   
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
   
   output$unalignedID <- renderText(as.character(
    datos$namedSequences[tabla_rows_selected(),]$ID))
   
   output$unalignedFASTA <- renderText(as.character(
    datos$namedSequences[tabla_rows_selected(),]$Sequence))
  }) # End bracket observeEvent entry selection