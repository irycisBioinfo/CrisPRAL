#----Running Pipeline----

observeEvent(input$Run_Example, {
 
 checks$example <- TRUE
 
 required_files <- c('../basic_Trials/Mutation_mus_R1_001.fastq', '../basic_Trials/Mutation_mus_R2_001.fastq', '../basic_trials/Reference.fasta',
                     '../basic_Trials/Primers.fasta', '../basic_Trials/Adapters-MiSeq.fasta')
 
 if(!any(file.exists(required_files))){
  
  return()
  
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
  
  system(paste('mkdir',datos$tmppipelinedir, sep = ' '))
  
  system(paste('echo ', R1_file(), '> /home/sergio/file.log', sep = ''))
  
  command = paste(
   CompletePATH,
   "/pipeline_v7.pl --r1 ",
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
   " --tmpdir ",
   as.character(datos$tmppipelinedir),
   collapse = "",
   sep = ""
  )
  
  system(command)#-Executes perl script
  
  incProgress( 1/5 ,detail = 'Alignments')
  
  if(!file.exists(paste(datos$tmppipelinedir,"/cluster.tsv", sep = ''))){
   #-Missing files
   
   output$Print2 <- renderText(':ERROR: Pipeline was broken; check files, extensions, reload and if problem persists contact administrator at sergio.fern1994@gmail.com')
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
  
  datos$allfasta <- read_delim(paste(datos$tmppipelinedir,"/final.fasta" ,sep = ''), 
                               delim = ' ', col_names =FALSE) %>% select(X1)
  
  datos$justfasta <- datos$allfasta %>% filter(!str_starts(X1, '>'))
  datos$justIds <- datos$allfasta %>% filter(str_starts(X1, '>'))
  datos$allfasta <- bind_cols(datos$justIds, datos$justfasta)
  colnames(datos$allfasta) <- c('ID','SEQ')
  datos$allfasta <- datos$allfasta %>% mutate_at(vars(ID), 
                                                 list(~if_else(str_starts(., '>'),
                                                               true = str_extract(., '(?<=>)[:graph:]+'), 
                                                               false = .)))
  
  datos$clusterF <- datos$cluster %>% select("ClusterN", "ID", "Identity")
  colnames(datos$clusterF) <- c('CLUSTER', 'ID', 'Identity')
  datos$clufast <- full_join(datos$allfasta, datos$clusterF) %>% group_by(CLUSTER) %>% nest() %>% arrange(CLUSTER)
  
  #datos$Tabla_raw is a rearrangement of datos$cluster to prepare the data for
  #fusing with Cluster-Sequences aswell as integrating with Alignment data
  
  datos$Tabla_raw = datos$cluster %>%
   group_by(ClusterN) %>%
   mutate(Abundance = n()) %>%
   ungroup() %>%
   mutate(Freq = 100 *Abundance / n()) %>%
   filter(Identity == "*") %>%
   select(ID, Abundance, Freq) %>%
   arrange(desc(Abundance))
  
  Total_Abundance = sum(datos$Tabla_raw['Abundance'])
  
  #----Associating clusterN with sequences----
  
  datos$Tabla_raw_clstrs <- datos$Tabla_raw %>% rowid_to_column('ClusterN') %>% 
   select(-Abundance, -Freq)
  
  datos$clustering <- full_join(datos$Tabla_raw_clstrs, datos$clufast %>% 
                                 unnest(), by = 'ID') %>% arrange(CLUSTER)
  datos$clustering <- datos$clustering %>% fill(ClusterN, .direction = 'down') %>% 
   select(-Identity, -CLUSTER ) %>% 
   group_by(ClusterN) %>% nest %>% 
   arrange(ClusterN)
  
  #Preparing for Alignment operations---- 
  #-asserting Reference does not provide an error:
  if(isFALSE(str_detect(Reference(), 'ERROR' ))){ #Checks whether Reference() returns an extension error.
   datos$Ref = readDNAStringSet(Reference())
  }else{
   if (!is.null(datos$Tabla)){
    Total_Abundance = sum(datos$Tabla_raw['Abundance'])
    output$Print2 <- renderText(paste("Total amount of Reads: ", Total_Abundance,                                           '
                                      :ERROR: File extension unknown; ', Reference(),'.
                                      Keeping previous results.', sep = ''))
    return()
    
   }else{
    
    output$Print2 <- renderText(Reference())
    return()}
   
  }
  
  datos$Sequences <- readDNAStringSet(paste(datos$tmppipelinedir,"/cluster", sep = ''))
  
  #Alignment for Sequences with Reference
  ls_ClusterAln_Ref = Clusters_Alignments(datos$Sequences, datos$Ref) #Function in Functions.R module.
  #Another Tabla variable is declared as an unsorted version together with
  #an associated Abundance for INDEL processing in graphing section.
  datos$Tabla_unsort <- ls_ClusterAln_Ref$tmp %>% rename(Deleted_bps = deletions) %>% rename(Inserted_bps = insertions)
  
  
  datos$aln = ls_ClusterAln_Ref$aln
  
  incProgress( 1/8 ,detail = 'Alignments')
  
  
  #Changing deletions and insertion columns to show total number of indels 
  #instead of total number of indel bps
  number_of_dels <- map(as.list(deletion(datos$aln)), ~length(.))
  total_deletions_per_cluster <- as_tibble(purrr::flatten_dbl(number_of_dels)) %>% rename(Deletions = value)
  number_of_ins <- map(as.list(insertion(datos$aln)), ~length(.))
  total_insertions_per_cluster <- as_tibble(purrr::flatten_dbl(number_of_ins)) %>% rename(Insertions = value)
  
  # datos$Tabla_unsort_total_indels <- cbind(cbind(datos$Tabla_unsort %>% select(-Deleted_bps, -Inserted_bps), total_deletions_per_cluster, total_insertions_per_cluster), datos$Tabla_unsort %>% select(Deleted_bps, Inserted_bps))
  
  datos$Tabla_unsort_total_indels <- cbind(datos$Tabla_unsort %>% select(-Deleted_bps, -Inserted_bps), total_deletions_per_cluster, datos$Tabla_unsort %>% select(Deleted_bps), total_insertions_per_cluster, datos$Tabla_unsort %>% select(Inserted_bps))
  
  datos$Tabla = inner_join(datos$Tabla_raw, datos$Tabla_unsort_total_indels) %>% 
   mutate(score = round(score,1), Freq = signif(Freq,2)) %>% 
   arrange(desc(Abundance))
  
  rownames(datos$Tabla) <- str_c('Cluster', rownames(datos$Tabla))
  
  datos$Tabla_Original <- datos$Tabla #datos$Tabla will be printed as a reactive function, line located before "Running Pipeline" block. datos$Tabla_Original is made to recover if datos$Tabla is altered by the user.
  output$tablaD <- renderDT(datos$Tabla, server = TRUE) #Table with multiple selection for PDF formatting
  # output$Print <- renderPrint(datos$Tabla[sort.default(rows_selected()),])
  output$Print <- renderDT(datos$Tabla[sort.default(rows_selected()),], selection = 'none')
  output$Print2 <- renderText(paste("Total amount of Reads: ", 
                                    Total_Abundance, sep = ''))
  
  #Data with named groups for FASTAS downloads.
  datos$namedSequences <- datos$Tabla %>% 
   left_join(datos$fasta) %>% select(ID, Sequence) %>%
   mutate(ID = paste(ID, '_', rownames_Tabla()[which(datos$Tabla == ID)] , sep = ''))
  datos$namedStringSet <- DNAStringSet(datos$namedSequences$Sequence)
  names(datos$namedStringSet) <- datos$namedSequences$ID
  
  #Is Target Introduced?
  if(!is.null(Target_FileInput()) && isFALSE(str_detect(Target(), 'ERROR' ))){ 
   incProgress( 1/8 ,detail = 'Looking for Target')
   datos$Target = readDNAStringSet(Target())
   
   #Alignment for Sequences with Reference
   ls_ClusterAln_Tar = Clusters_Alignments(datos$Sequences, datos$Target) #Function in Functions.R module.
   #Another Tabla variable is declared as an unsorted version together with
   #an associated Abundance for INDEL processing in graphing section.
   datos$TablaT_unsort <- ls_ClusterAln_Tar$tmp %>% rename(Deleted_bps = deletions) %>% rename(Inserted_bps = insertions)
   datos$unsortT_ID_Abundance <- left_join(datos$TablaT_unsort %>% select(ID),
                                           datos$Tabla_raw %>% select(-Freq),  by = 'ID') #What is this used for?
   
   incProgress( 1/8 ,detail = 'Alignments... renaming')
   
   datos$alnT = ls_ClusterAln_Tar$aln
   
   incProgress( 1/8 ,detail = 'Alignments')
   
   
   #Changing deletions and insertion columns to show total number of indels 
   #instead of total number of indel bps
   number_of_delsT <- map(as.list(deletion(datos$alnT)), ~length(.))
   total_deletions_per_clusterT <- as_tibble(purrr::flatten_dbl(number_of_delsT)) %>% rename(Deletions = value)
   
   number_of_insT <- map(as.list(insertion(datos$alnT)), ~length(.))
   total_insertions_per_clusterT <- as_tibble(purrr::flatten_dbl(number_of_insT)) %>% rename(Insertions = value)
   
   datos$TablaT_unsort_total_indels <- cbind(cbind(datos$TablaT_unsort %>% 
                                                    select(-Deleted_bps, -Inserted_bps), total_deletions_per_clusterT, total_insertions_per_clusterT), 
                                             datos$TablaT_unsort %>% 
                                              select(Deleted_bps, Inserted_bps))
   
   datos$Tabla_Target = inner_join(datos$Tabla_raw, datos$TablaT_unsort_total_indels) %>% 
    mutate(score = round(score,1), Freq = signif(Freq,2)) %>% 
    arrange(desc(Abundance))
   
   output$Target_Location <- renderText(paste(Target_location()[2], Target_location()[1], sep = ' '))
   
   rownames(datos$Tabla_Target) <- str_c('Cluster', rownames(datos$Tabla_Target))
   
   datos$Tabla_Target_Original <- datos$Tabla_Target
   # output$tablaT <- renderDT(datos$Tabla_Target, editable = TRUE, selection = "single", 
   #                           options = list(search = list(smart = FALSE)), server = TRUE)
   
  }else if(is.null(Target_FileInput())){
   output$Target_Location <- renderText(paste('No Target was introduced', 
                                              sep = ' '))
   
  }else if(str_detect(Target(), 'ERROR' )){ #Checks whether Target() return an extension error.
   
   output$Target_Location <- renderText(Target())
  }
 })
})