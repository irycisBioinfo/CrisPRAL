# Server-Mosaic

########################### DESCRIPTION ##################################
# Loads script for pipeline execution with the stablished parameters and
# loads the main Tables of results. Also defines code to allow user edits on table.
##########################################################################.

#----Running Pipeline----

pipeline <- reactive({if(input$reverse_complement){
  return("/Scripts/pipeline_vUMI.pl --r1 ")
  }else{return("/Scripts/pipeline_v7.pl --r1 ")}
  })

observeEvent(input$Accept, {
 checks$example <- FALSE
 
 required_values <- list("Reference sequence" = input$Reference$name, 'Fastq File' = input$R1$name)
 name <- names(required_values[required_values == 'NULL'])
 
 if( any(required_values == 'NULL') ) {
  
  name <- names(required_values[required_values == 'NULL'])
  
  showModal(modalDialog(
   fluidPage(
    h1(strong("Warning"),align="center"),
    hr(),
    fluidRow(column(4,offset = 4,
                    div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                        img(src="caution-icon.png", height='150', width='150', align="middle")))),
    h2(paste("A ", name ," is required to align reads." ,sep = ''), align = 'center'),
    easyClose = TRUE,
    footer = NULL
   )))
  
 }
 
 req(input$R1, input$Reference)
 if(input$single_end == FALSE){
  
  if(Check_pair_input_files() == FALSE){
   
   return()
   
  }
 }
 
 #-Activates Clustering and Alignment tab
 updateTabsetPanel(session, 'Visualization', selected = 'Clusterization_&_Alignment')
 
 datos$Tabla_prev <- NULL
 datos$Pie_data <- NULL
 output$Pie_summary <- renderPlotly({print(Edit_Pie_chart())})
 output$Mut_Freq <- renderPlotly({print(NULL)})
 output$Del_sizes <- renderPlotly({print(NULL)})
 output$In_sizes <- renderPlotly({print(NULL)})
 output$Del_loc <- renderPlotly({print(NULL)})
 output$In_loc <- renderPlotly({print(NULL)})
 
 datos$text <- 'Click on a Cluster to show its alignment.'
 output$Display_Unaligned <- renderUI({})
 output$Blast_Search <- renderUI({})
 output$downloadFASTA_file <- renderUI({})
 output$downloadPAIR_file <- renderUI({})
 output$Unaligned_Fasta <- renderUI({})
 output$score = renderText(paste(''))
 
 if(file.exists('default_charts.Rdata')){
  
  system(paste('rm default_charts.Rdata', sep = ''))
  system(paste('rm default_charts_doc.Rdata', sep = ''))
  
 }
 
 
 withProgress(message = "Performing analisis: ", detail = 'Clustering, may take a while...', {
  
  datos$tmppipelinedir <- paste(session$token, '/data', str_sub(tempfile(), start = 21) ,sep = '')
  
  dir.create(datos$tmppipelinedir)
  
  command = paste(
   CompletePATH,
   as.character(pipeline()),
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
   " --primer-error-rate ",
   input$primer_error_rate,
   " --tmpdir ",
   as.character(datos$tmppipelinedir),
   collapse = "",
   sep = ""
  )
  
  system(command) #-Executes perl script
  
  incProgress( 1/5 ,detail = 'Alignments')
  
  if(!file.exists(paste(datos$tmppipelinedir,"/cluster.tsv", sep = ''))){
   #-Missing files
   
   output$Print2 <- renderText(paste0(':ERROR: Pipeline was broken; check files, extensions, reload and if problem persists contact administrator at sergio.fern1994@gmail.com'))
   return()
   
  }
  
  #This directory has to be reactive
  datos$fasta <- read_tsv(paste(datos$tmppipelinedir,"/cluster.tsv", sep = ''), 
                          col_names = FALSE)
  datos$cluster <- read_tsv(paste(datos$tmppipelinedir,"/cluster.bak.tsv", sep = ''), 
                            col_names = FALSE)
  
  if (isEmpty(colnames(datos$fasta))){#-If dataset is completely filtered 
   #due to any restriction app crashes
   if (!is.null(datos$Tabla)){
    Total_Abundance = sum(datos$Tabla_raw['Abundance'])
    output$Print2 <- renderText(paste("Total amount of Reads: ", Total_Abundance, 
                                      '
                                      :WARNING: All data of new process was filtered. Posibly something happened, adjust parameters & make sure files are correctly introduced.
                                      Keeping previous Results...', sep = ''))
    return()
    
   }else{
    output$Print2 <- renderText(':WARNING: All data was filtered. Posibly something happened, adjust parameters & make sure files are correctly introduced.')
    return()
   }
   
  }else{
   colnames(datos$fasta) <- c("ID", "Sequence")
   colnames(datos$cluster) <- c("ClusterN", "Length", "ID", "Identity")
  } 

  allfasta <- read_delim(paste(datos$tmppipelinedir,"/final.fasta" ,sep = ''), 
                               delim = ' ', col_names =FALSE) %>% select(X1)
  
  justfasta <- allfasta %>% filter(!str_starts(X1, '>'))
  justIds <- allfasta %>% filter(str_starts(X1, '>'))
  allfasta <- bind_cols(justIds, justfasta)
  colnames(allfasta) <- c('ID','SEQ')
  allfasta <- allfasta %>% mutate_at(vars(ID), 
                                                 list(~if_else(str_starts(., '>'),
                                                               true = str_extract(., '(?<=>)[:graph:]+'), 
                                                               false = .)))
  
  # We want to double the amount of total reads when coming from a paired-end experiment.
  if(input$single){ correct_paired_reads = 1 }else{ correct_paired_reads = 2 }
  
  
  datos$Tabla_raw = datos$cluster %>%
   group_by(ClusterN) %>%
   mutate(Abundance = n()*correct_paired_reads) %>%
   ungroup() %>%
   mutate(Freq = 100 *Abundance / (n()*correct_paired_reads)) %>%
   filter(Identity == "*") %>%
   select(ID, Abundance, Freq, Length) %>%
   arrange(desc(Abundance))
  
  Total_Abundance = sum(datos$Tabla_raw['Abundance'])
  
  #----Associating clusterN with sequences----
  
  # datos$Tabla_raw_clstrs <- datos$Tabla_raw %>% rowid_to_column('ClusterN') %>% 
  #  select(-Abundance, -Freq)
  
  # datos$clustering <- full_join(datos$Tabla_raw_clstrs, datos$clufast %>% unnest(), by = 'ID') %>% arrange(CLUSTER)
  # datos$clustering <- datos$clustering %>% fill(ClusterN, .direction = 'down') %>% 
  #  select(-Identity, -CLUSTER ) %>% 
  #  group_by(ClusterN) %>% nest %>% 
  #  arrange(ClusterN)
  
  #Preparing for Alignment operations---- 
  #-asserting Reference does not provide an error:
  if(isFALSE(str_detect(Reference(), 'ERROR' ))){ #Checks whether Reference() returns an extension error.
   datos$Ref = readDNAStringSet(Reference())
  }else{
   if (!is.null(datos$Tabla)){
    Total_Abundance = sum(datos$Tabla_raw['Abundance'])
    output$Print2 <- renderText(paste("Total amount of Reads: ", Total_Abundance, '
                                      :ERROR: File extension unknown; ', Reference(),'.
                                      Keeping previous results.', sep = ''))
    return()
    
   }else{
    
    output$Print2 <- renderText(Reference())
    return()}
   
  }
  
  datos$Sequences <- readDNAStringSet(paste(datos$tmppipelinedir,"/cluster", sep = ''))
  #Alignment for Sequences with Reference
  ls_ClusterAln_Ref = Clusters_Alignments(datos$Sequences, datos$Ref, 
                                          gapOpening = input$gap_open, 
                                          gapExtension = input$gap_extend,
                                          type = input$general_allType) #Function in Functions.R module.
  #Another Tabla variable is declared as an unsorted version together with
  #an associated Abundance for INDEL processing in graphing section.
  datos$Tabla_unsort <- ls_ClusterAln_Ref$tmp
  
  
  ############.
  # alignment data is extracted and stored for the whole app since serverGraphs.R uses this to produce the graphic representation
  datos$aln <- ls_ClusterAln_Ref$aln
  ############.
  
  incProgress( 1/8 ,detail = 'Alignments')
  
  #- Finetune datos$Tabla
  datos$Tabla = inner_join(datos$Tabla_raw, ls_ClusterAln_Ref$tmp) %>% 
   mutate(score = round(score,1), Freq = signif(Freq,2)) %>% 
   arrange(desc(Abundance))# %>% select(-width,-start,-end)
  
  ###########################################################.
  ##### Mismatch extraction and joining for datos$Tabla #####
  ###########################################################.
  
  datos$Tabla <- mismatches_break_down(ls_ClusterAln_Ref, datos$Tabla)
  
  ##################################.
  #### OUTPUTING Reference ####
  ##################################.
  
  rownames(datos$Tabla) <- str_c('Cluster', rownames(datos$Tabla))
  
  datos$Tabla_Original <- datos$Tabla #datos$Tabla will be printed as a reactive function, line located before "Running Pipeline" block. datos$Tabla_Original is made as recovery if datos$Tabla is altered by the user.
  output$tablaD <- renderDT(Tabla(),filter = 'top', server = TRUE) #Table with multiple selection for PDF formatting
  output$PrintD <- renderDT(Tabla()[sort.default(rows_selected()),], selection = 'none')
  output$Print2 <- renderText(paste("Total amount of Reads: ", 
                                    Total_Abundance, sep = ''))
  
  ###########################################################.
  ##### FASTA data preparation #####
  ###########################################################.
  
  datos$namedSequences <- datos$Tabla %>% 
   left_join(datos$fasta) %>% select(ID, Sequence) %>%
   mutate(ID = paste(ID, '_', rownames_Tabla()[which(datos$Tabla == ID)] , sep = ''))
  datos$namedStringSet <- DNAStringSet(datos$namedSequences$Sequence)
  names(datos$namedStringSet) <- datos$namedSequences$ID
  
#########################################################################.
##### Is Target Introduced? #####
#########################################################################.
  
  if(!is.null(Target_FileInput()) && isFALSE(str_detect(Target(), 'ERROR' ))){ 
   incProgress( 1/8 ,detail = 'Looking for Target')
   datos$Target = readDNAStringSet(Target())
   
   #Alignment for Sequences with Target
   ls_ClusterAln_Tar = Clusters_Alignments(datos$Sequences, datos$Target,
                                           gapOpening = input$gap_open, 
                                           gapExtension = input$gap_extend,
                                           type = input$general_allType) #Function in Functions.R module.
   
   datos$Tabla_Target = inner_join(datos$Tabla_raw, ls_ClusterAln_Tar$tmp) %>% 
    mutate(score = round(score,1), Freq = signif(Freq,2)) %>% 
    arrange(desc(Abundance))# %>% select(-width,-start,-end) # Remove columns deemed too confusing for a user.
   
   ##################################################################.
   ##### Mismatch extraction and joining for datos$Tabla_Target #####
   ##################################################################.
   
   datos$Tabla_Target <- mismatches_break_down(ls_ClusterAln_Tar, datos$Tabla_Target)
   
   output$Target_Location <- renderText(paste(Target_location()[1], Target_location()[2], sep = ' '))
   
   rownames(datos$Tabla_Target) <- str_c('Cluster', rownames(datos$Tabla_Target))
   
   datos$Tabla_Target_Original <- datos$Tabla_Target
   
   ##########################.
   #### OUTPUTING Target ####
   ##########################.
   
   # Declaring Target Table for Download tab.
   output$tablaTD <- renderDT(TablaT(), server = TRUE) #Table with multiple selection for PDF formatting
   output$PrintTD <- renderDT(TablaT()[sort.default(rows_selected()),], selection = 'none')
   
  }else if(is.null(Target_FileInput())){
   output$Target_Location <- renderText(paste('No Target was introduced', 
                                              sep = ' '))
   
  }else if(str_detect(Target(), 'ERROR' )){ #Checks whether Target() return an extension error.
   
   output$Target_Location <- renderText(Target())
  }
 })
}) # observeEvent close bracket

output$temp_text <- renderUI({
 verbatimTextOutput('explicative_text')
})

# Rownames can be changed and these affect the outputing of the pie chart.
rownames_Tabla <- reactive({
 rownames(datos$Tabla)
}) 
###############################################################################.
##### User Table Edit ####
###############################################################################.

  ###########################.
  #-Table Filtering ####
  ###########################.

observeEvent(input$Filter, {
 #-Table Reference:
 datos$Tabla <- datos$Tabla_Original %>% rownames_to_column(var = 'rowname') %>% 
  filter(score >= input$score_threshold) %>% 
  filter(Abundance >= input$abundance_minimum) %>%
  mutate_at(vars(Freq), funs(signif(100*Abundance/ sum(Abundance),2))) %>% 
  column_to_rownames('rowname')
 
 #-Table Target:
 if(!is.null(datos$Tabla_Target_Original)){
  datos$Tabla_Target <- datos$Tabla_Target_Original %>% rownames_to_column(var = 'rowname') %>% 
   filter(score >= input$score_threshold) %>% 
   filter(Abundance >= input$abundance_minimum) %>%
   mutate_at(vars(Freq), funs(signif(100*Abundance/sum(Abundance),2))) %>% 
   column_to_rownames('rowname')
  
 }
})
#-Back to Default values:
observeEvent(input$Default, {
 datos$Tabla <- datos$Tabla_Original
 if(!is.null(datos$Tabla_Target_Original)){
  datos$Tabla_Target <- datos$Tabla_Target_Original}
}
)

  ################################.
  #-Rowname Customization ####
  ################################.

proxy_tablaR = dataTableProxy('tablaR')
proxy_tablaT = dataTableProxy('tablaT')
proxy_tablaD = dataTableProxy('tablaD')

observeEvent(c(input$tablaR_cell_edit, input$tablaT_cell_edit),{
 
 edit <- list(tablaR = input$tablaR_cell_edit, tablaT = input$tablaT_cell_edit)
 
 info = edit[edit != 'NULL']
 
 # str(info) Check info for debugging porpouses
 
 # Using switch with input$Alignment to determine which table has been edited. 
 # Necessary since input$*_cell_edit is read_only and lingers
 
 Alignment <- switch(input$Alignment, 'Reference' = 'tablaR', 'Target' = 'tablaT')
 
 i = info[[Alignment]]$row 
 j = info[[Alignment]]$col 
 v = info[[Alignment]]$value
 
 if(j == 0){
  rownames(datos$Tabla)[i] <- v
  replaceData(proxy_tablaR, datos$Tabla, resetPaging = FALSE) # important
  replaceData(proxy_tablaD, datos$Tabla, resetPaging = FALSE)
  
  idx_tranform_to_originalR <- grep(datos$Tabla_Original$ID, pattern = datos$Tabla[1,i])
  
  rownames(datos$Tabla_Original)[idx_tranform_to_originalR] <- v
  
  if(!is.null(datos$Tabla_Target)){
   rownames(datos$Tabla_Target)[i] <- v
   replaceData(proxy_tablaT, datos$Tabla_Target, resetPaging = FALSE) # important
   
   idx_tranform_to_originalT <- grep(datos$Tabla_Target_Original$ID, pattern = datos$Tabla_Target[1,i])
   rownames(datos$Tabla_Target_Original)[idx_tranform_to_originalT] <- v
   
  }
  
  if(!is.null(charts$PC_online)){
   output$Pie_summary <- renderPlotly({print(Edit_Pie_chart())})
  }
  
 }else{
  datos$Tabla[i, j] <<- DT::coerceValue(datos$Tabla[i, j], datos$Tabla[i, j])
  replaceData(proxy_tablaR, datos$Tabla, resetPaging = FALSE) # important
  
  showModal(modalDialog(
   
   fluidPage(
    h1(strong("Warning"),align="center"),
    hr(),
    fluidRow(column(4,offset = 4,
                    div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                        img(src="caution-icon.png", height='100', width='100', align="middle")))),
    h2("Only first column is editable", align = 'center'),
    easyClose = TRUE,
    footer = NULL
   )
   
  ))
 }
 reset('tablaR_cell_edit')
 reset('tablaT_cell_edit')
})