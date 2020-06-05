# serverMosaic <- function(input, output, session) {
 #---Mosaic Analyser SERVER----
  #-Reset block-
  
   reset('R1')
   reset('R2')
   reset('Reference')
   reset('Target')
   
   observeEvent(input$toggle2, {
    toggleDropdownButton(inputId = "mydropdown")
     }, ignoreInit = TRUE)
   
   #------TESTING-----
   
   # output$coloured_text <- renderUI({
   #  wellPanel(
   #   tags$div(
   #    HTML(paste("This text is ", 
   #             tags$span(style="background-color:rgba(255,0,0,0.40)", "r"),
   #             tags$span(style="background-color:rgba(0,255,0,0.40)", "e"),
   #             tags$span(style="background-color:rgba(0,0,255,0.40)", "d"), sep = ""))
   # ))})
   
  #-Compartimentalization line-
  system(paste('mkdir', session$token, sep = ' '))

  rows_selected <- reactive({ if( length(input$tablaD_rows_selected ) != 0 ){ 
   return(input$tablaD_rows_selected )
  }else if ( input$all_clusters ){ 
   return(c(1:length(datos$Sequences)))
  }else{ return(c(1:10)) }})
  
  warningCLSTR <- reactive({if(length(datos$clustersStringSet) > 100){
     return('Warning: Number of sequences for Multiple Sequence Alignment exceed max computing limit, 100 will be computed but only the first 50 entries will be shown')
    }else if(length(datos$clustersStringSet) > 50){
     return('Warning: Number of sequences for Multiple Sequence Alignment exceed max printing limit, only the first 50 entries will be shown')
    }else{return(FALSE)}
   })
  
  #---Reactive values-Data storage----  
  
  datos <- reactiveValues()
  charts <- reactiveValues()
  file_upload <- reactiveValues(upload_Adapter_state = NULL)
  file_upload <- reactiveValues(upload_Primer_state = NULL)
  
  #-Reactive input to allow for compressed files processing.
  R1_file <- reactive({
    if (str_sub(input$R1$datapath, -3) == '.gz'){
      #Function "Unzip_file" defined in Functions.R
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
      extensionR <- str_split(input$Reference$datapath, '\\.')[[1]][length(str_split(input$Reference$datapath, '\\.')[[1]])]
      lastIndexR <- nchar(extensionR)
      T_Refence_dir = paste(str_sub(input$Reference$datapath, 1,-(lastIndexR+3)),
                            'Filtered_Reference.fasta', sep = '')
      
      if(extensionR == 'txt'){
       RefAsFasta <- paste(str_sub(input$Reference$datapath, 1,-(lastIndexR+3)) ,'RefAsFasta.fasta', sep = '')
       system(paste( 'cp', input$Reference$datapath, RefAsFasta, sep = ' '))
      }else if(extensionR == 'fasta'){
       RefAsFasta <- input$Reference$datapath
      }else{
       return(':ERROR: Unknown extension, please make sure that the Reference file has extension .txt or .fasta')
      }
      
      system(paste('cutadapt -n 2 -g ', 
                   as.character(PrimerR1_5()), ' -a ', 
                   as.character(PrimerR1_3()), ' -o ', 
                   T_Refence_dir,' ',
                   RefAsFasta, sep = '' ))
      return(T_Refence_dir)
    }
    else{
      return(input$Reference$datapath)
    }
  })
  
  Target_FileInput <- callModule(inputReset_server, 'reset_Target')
  Target <- reactive({
    if(file_upload$upload_Primer_state == 'uploaded' && !is.null(Target_FileInput())){
     #This following code is necessary to process both .fasta & .txt extensions.
     extensionT <- str_split(Target_FileInput()$datapath, '\\.')[[1]][length(str_split(Target_FileInput()$datapath, '\\.')[[1]])]
     lastIndexT <- nchar(extensionT)
      T_Target_dir = paste(str_sub(Target_FileInput()$datapath, 1,-(lastIndexT+3)),
                           'Filtered_Target.fasta', sep = '')
      
      if(extensionT == 'txt'){
       TarAsFasta <- paste(str_sub(Target_FileInput()$datapath, 1,-(lastIndexT+3)) ,'TarAsFasta.fasta', sep = '')
       system(paste( 'cp', Target_FileInput()$datapath, TarAsFasta, sep = ' '))
      }else if(extensionT == 'fasta'){
       TarAsFasta <- Target_FileInput()$datapath
      }else{
       return(':ERROR: Unknown extension; please make sure the Target file has a .txt or .fasta')
      }
      
      system(paste('cutadapt -n 2 -g ', 
                   as.character(PrimerR1_5()), ' -a ', 
                   as.character(PrimerR1_3()), ' -o ', 
                   T_Target_dir, ' ', 
                   TarAsFasta, sep = '' ))
      
      return(T_Target_dir)
    }else{
      return(Target_FileInput()$datapath)
    }
  })
  
  ObsT <- observe({
    
    if(!is.null(Target_FileInput())){
      showTab( inputId ='Alignment', target = 'Target' )
    }else{
      hideTab( inputId ='Alignment', target = 'Target')
    }})
  
  #-Table data made reactive.
  Tabla <- reactive({datos$Tabla})
  TablaT <- reactive({datos$Tabla_Target})
  
  Target_location <- reactive({find_target_location(Tabla(), TablaT())})
  
  #Generates reactive UI depending on the presence of a Target or not.
  output$is.nullTarget1 <- renderUI({
    selectInput(
      "alnTo",
      "Alignment to:",
      if (!is.null(Target_FileInput())){c("Reference","Target")
        }else{c("Reference")})
  })
  
  #----Filtering Parameters declarations----
  
  #This following set of observeEvents allows the user to change his mind on 
  #uploaded Adapter files.
  observeEvent(c(input$AdapterFile, input$uploadAdapters), {
    file_upload$upload_Adapter_state <- 'uploaded'
  })
  observeEvent(c(input$resetA, input$input_type_A), {
    file_upload$upload_Adapter_state <- 'reset'
  })
  
  Adapter_Directinput <- eventReactive(input$uploadAdapters,
                                       { c(input$R1_Adapter,input$R2_Adapter) })
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
      if(input$input_type_A == 'File_input_A'){
        return(readDNAStringSet(Adapter_FileInput()$datapath))
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
  
  output$checkR1A <- renderText(paste("R1 Adapter: ",AdapterR1(), sep = ''))
  output$checkR2A <- renderText(paste("R2 Adapter: ",AdapterR2(), sep = ''))
  
  #Same than the Adapter strings but with the Primers.
  observeEvent(c(input$PrimerFile, input$uploadPrimers), {
    file_upload$upload_Primer_state <- 'uploaded'
  }) 
  observeEvent(c(input$resetP, input$input_type_P), {
    file_upload$upload_Primer_state <- 'reset'
  })
  
  Primer_Directinput <- eventReactive(input$uploadPrimers,{ 
    c(input$F_Primer1,input$R_Primer2) 
    })
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
  #Primers loaded as reactive functions to prevent residual data lingering.
  Primer5_sequences <- reactive({
    
    if (is.null(Primer_FileInput())){
      return("Empty")
    }else if (!is.null(Primer_FileInput())){
      if(input$input_type_P == 'File_input_P'){return(readDNAStringSet(
                                                      Primer_FileInput()$datapath))
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
  
  output$Explicative_Text1 <- renderText(paste(
    'The file should contain four lines, two', 
    'lines with ">" indicating Primer ID and', 
    'another two lines with the sequences', 'themselves. Ex:',
    '>Primer1', 'ACGT', '>Primer2', 'GTAC',sep = '\n'))
  
  output$checkFP1 <- renderText(paste("Read 1, 5' end: ",PrimerR1_5(), 
                                      sep = ''))
  output$checkFP2 <- renderText(paste("Read 1, 3' end: ",PrimerR1_3(), 
                                      sep = ''))
  output$checkRP1 <- renderText(paste("Read 2, 5' end: ",PrimerR2_5(), 
                                      sep = ''))
  output$checkRP2 <- renderText(paste("Read 2, 3' end: ",PrimerR2_3(), 
                                      sep = ''))
  #---------Modal dialog WARNING----
  
  Check_pair_input_files <- reactive({ #Allow warning to be disabled
   
   if(is.null(input$R2)){
    
    showModal(modalDialog(
     fluidPage(
      h1(strong("Warning"),align="center"),
      hr(),
      fluidRow(column(4,offset = 4,
                      div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                          img(src="caution-icon.png", height='150', width='150', align="middle")))),
      h2(paste("Two files are required to align reads in paired-end mode." ,sep = ''), align = 'center'),
      easyClose = TRUE,
      footer = NULL
     )))
    
   }
   
   req(input$R2)
   req(input$R1)
   
   R1 <- substring(input$R1$name, seq(1, nchar(input$R1$name), 1), seq(1, nchar(input$R1$name), 1))
   R2 <- substring(input$R2$name, seq(1, nchar(input$R2$name), 1), seq(1, nchar(input$R2$name), 1))
   
   if( sum(R1 != R2) > 1){
    showModal(modalDialog(
     fluidPage(
      h1(strong("Warning"),align="center"),
      hr(),
      fluidRow(column(4,offset = 4,
                      div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                          img(src="caution-icon.png", height='150', width='150', align="middle")))),
      h2("File names do not seem to match a paired-end experiment", align = 'center'),
      hr(),
      h4("Please make sure files are correct before proceeding.", align = 'center'),
      easyClose = TRUE,
      footer = NULL
     )))
    
    return(FALSE)
    
   }else if(sum(R1 != R2) == 0){
    
    showModal(modalDialog(fluidPage(
     h1(strong("Warning"),align="center"),
     hr(),
     fluidRow(column(4,offset = 4,
                     div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                         img(src="caution-icon.png", height='150', width='150', align="middle")))),
     h2("File names do not seem to match a paired-end experiment", align = 'center'),
     hr(),
     h4("Please make sure files are correct before proceeding.", align = 'center'),
     easyClose = TRUE,
     footer = NULL
    )))
    
    return(FALSE)}
   
   else{return(TRUE)}
   
  })
  
  #----Reactive table printing----
  
  output$tablaR <- renderDT(Tabla(), selection = "single", editable = TRUE, 
                            options = list(search = list(smart = FALSE)), server = TRUE)
  
  output$tablaT <- renderDT(TablaT(), editable = TRUE, selection = "single", 
                            options = list(search = list(smart = FALSE)), server = TRUE)
  
  #----Running Pipeline----
  
  observeEvent(input$Accept, {
   
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
    
    updateTabsetPanel(session, 'Visualization', selected = 'Clusterization_&_Alignment')
    
    datos$Pie_data <- NULL
    output$Pie_summary <- renderPlotly({print(Edit_Pie_chart())})
    output$Mut_Freq <- renderPlotly({print(NULL)})
    output$Del_sizes <- renderPlotly({print(NULL)})
    output$In_sizes <- renderPlotly({print(NULL)})
    output$Del_loc <- renderPlotly({print(NULL)})
    output$In_loc <- renderPlotly({print(NULL)})
    
    withProgress(message = "Performing analisis: ",
                 detail = 'Clustering, may take a while...', {
      #-Assistive reactive functions
      datos$tmppipelinedir <- paste(session$token, '/data', str_sub(tempfile(), start = 21) ,sep = '')
      
      system(paste('mkdir',datos$tmppipelinedir, sep = ' '))
      
      command = paste(
        CompletePATH,
        "/pipeline_v6_min_qual_pipedir.pl --r1 ",
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
      #-asserting Reference does not provide an error.
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
      datos$unsort_ID_Abundance <- left_join(datos$Tabla_unsort %>% select(ID), 
                                             datos$Tabla_raw %>% select(-Freq),  by = 'ID')
      
      incProgress( 1/8 ,detail = 'Alignments... renaming')
      
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
        
        output$Target_Location <- renderText(paste(Target_location()[2], paste(Target_location()[1], collapse = ', '), sep = ' '))
        
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
  
  output$temp_text <- renderUI({
    verbatimTextOutput('explicative_text')
   })
  
  #-Alignment and Table filtering dinamic UI-
  output$alignment_with_filtering <- renderUI({
   
   if(is.null(datos$Tabla_Original)){
    
    return(output$explicative_text <- NULL)
    
   }
   
   fluidRow(
    column(5,
           wellPanel(
            h3("Alignment"),
            selectInput(
             "alnType",
             "Alignment method",
             c("global", "local", "overlap", "global-local", "local-global")
            ),
            uiOutput('is.nullTarget1'),
            verbatimTextOutput('Print2'),
            uiOutput('coloured_text'),
            verbatimTextOutput('Target_Location'))),
    column(5, 
           wellPanel(h3("Cluster filtering"),
                     sliderInput('score_threshold',label = p('Scoring threshold'), 
                                 min = round(min(datos$Tabla_Original %>% select(score))),
                                 max = round(max(datos$Tabla %>% select(score))), 
                                 value = round(min(datos$Tabla %>% select(score)))),
                     fluidRow(
                      column(2,
                             actionButton('Filter', 'Accept')),
                      column(2, actionButton('Default', 'Default'))))
    )
   )})

  
 rownames_Tabla <- reactive({
  rownames(datos$Tabla)
 }) 
   
 #----User Table Edit----
  #-Table Filtering:
 observeEvent(input$Filter, {
  #-Table Reference:
  datos$Tabla <- datos$Tabla_Original %>% rownames_to_column(var = 'rowname') %>% 
   filter(score >= input$score_threshold) %>% 
   mutate_at(vars(Freq), funs( round(100*Abundance/ sum(Abundance)))) %>% 
   column_to_rownames('rowname')

  #-Table Target:
  if(!is.null(datos$Tabla_Target_Original)){
  datos$Tabla_Target <- datos$Tabla_Target_Original %>% rownames_to_column(var = 'rowname') %>% 
   filter(score >= input$score_threshold) %>% mutate_at(vars(Freq), funs(round(100*Abundance/sum(Abundance)))) %>% 
   column_to_rownames('rowname')

  }
 })
   #-Back to Default values:
 observeEvent(input$Default, {
  datos$Tabla <- datos$Tabla_Original
  if(!is.null(datos$Tabla_Target_Original)){
  datos$Tabla_Target <- datos$Tabla_Target_Original}})

 #-Rowname Custimization:
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
observeEvent( c( input$tablaR_rows_selected, input$tablaT_rows_selected, 
               input$alnType, input$alnTo), {
                
     req( input$tablaR_rows_selected != datos$prev_selected_row || input$alnType != datos$prev_alnType || input$alnTo != datos$prev_alnTo)
                 
   withProgress( message = "Performing analisis", {
     
     req( datos$Tabla )
     req( tabla_rows_selected( ) )

     
     datos$prev_selected_row <- tabla_rows_selected()
     
     Seq <- datos$namedStringSet[tabla_rows_selected()]
     
     if(input$alnTo == "Reference")
     {
       datos$aln2 = pairwiseAlignment( datos$Ref, Seq, 
                                      type = input$alnType, gapOpening=10, 
                                      gapExtension=5 )
     }else{
       datos$aln2 = pairwiseAlignment( datos$Target, Seq, 
                                      type = input$alnType, gapOpening=10, 
                                      gapExtension=5 )
     }
     
     comp = character(length( datos$aln2@pattern ))
     ref = unlist(strsplit(as.character( datos$aln2@pattern ), 
                           split = "" ))
     query = unlist(strsplit(as.character( datos$aln2@subject ), 
                           split = "" ))
     
     for (i in 1:length( ref ))
     {
       if (ref[i] != query[i])
       {
         comp[i] = query[i]
       } else {
         comp[i] = "."
       }
     }
     
     pos <- Position_vector(as.character(datos$aln2@pattern), 
                            ref) #Generates positioning vector
     
     if (!is.null(tabla_rows_selected()))
     {
       datos$text = as.data.frame(c(
         #Not working well enough, fails sometimes
         paste(pos, sep = '', collapse = ''), 
         paste(ref, sep = "", collapse = ""),
         paste(query, sep = "", collapse = ""),
         paste(comp, sep = "", collapse = "")
       ))
       
       if(input$alnTo == "Reference")
       {
         rownames(datos$text) = c("Position", "Reference", 
                          datos$Tabla$ID[tabla_rows_selected()], 
                          "Comparisson")
       }else{
         rownames(datos$text) = c("Position", "Target", 
                          datos$Tabla$ID[tabla_rows_selected()], 
                          "Comparisson")
       }
       
       colnames(datos$text) = c("Alignment")
       output$alignment = renderPrint(datos$text)
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

  #-Reactive functions for CallModule and downloads----
  
  aln2PAIR <- reactive({datos$aln2})
  aln2FASTA <- reactive({alignedSubject(datos$aln2)})
  
  namedStringSet <- reactive({
   names(datos$namedStringSet) <- str_c(str_split(names(datos$namedStringSet), '_')[[1]][1], 
                                        paste('_', rownames_Tabla(), sep = ''))
   return(datos$namedStringSet)})
  
  clustersStringSet <- reactive({datos$clustersStringSet})
  clustersMSA <- reactive({
   if (length(datos$clustersStringSet) <= 100){
    msa(datos$clustersStringSet)
    }else if(length(datos$clustersStringSet) >= 100){
     msa(datos$clustersStringSet[1:100])
    }
   
   })
  
  
  prettyprintCrashpreventor <- reactive({
     if (nrow(clustersMSA())/50 < 1){
      
     c(1:nrow(clustersMSA()))
      
     }else if(nrow(clustersMSA())/50 > 1){
      
     unique(round(seq.int(from = 1, to = nrow(clustersMSA()), length.out = 50)))
      
     }
  })
  
  input_filename_generic <- reactive({paste(
                                      str_replace(#Removes posible double 
                                                  #underscore "_"
                                      str_replace(#Removes Read number ID
                                        substring(text = input$R1$name, 
                                                  first = 1, 
                                                  last = nchar(input$R1$name)-(sum(nchar(str_split(input$R1$name, '\\.')[[1]][-1]))+length(str_split(input$R1$name, '\\.')[[1]][-1])))
                                        ,'_R1',''), '__',''), sep='')})
  input_filename_results <- reactive({paste( input_filename_generic(), 
                                             '_Results', sep='' )})
  input_filename_group <- reactive({paste( input_filename_generic(), 
                                           '_Group', tabla_rows_selected(), 
                                           sep='')})
  
  input_filename_pair <- reactive({paste( input_filename_group() , 
                                          '_pairwise_alignment', sep='')})
  
  input_filename_clstrs <- reactive({paste(input_filename_group(), 
                                           '_clstrs' ,sep = '')})
  input_filename_clstrs_msa <- reactive({paste(input_filename_clstrs(), '_msa', sep = '')})
  
  
  
  datos$alignmentdir <- 'empty' #Temporal asignment
   
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste(input_filename_results(),'pdf', sep = '.')
    },
    content = function(file) {
      src <- normalizePath('reports.Rmd')
      out <- rmarkdown::render( 'reports.Rmd', 
                                params = list(selection = rows_selected(),
                                             table = datos$Tabla %>% rownames_to_column() %>% 
                                              rename(Tag = 'rowname'),
                                             chartPie = charts$tmpFilePiepng,
                                             chartDS = charts$tmpFileDSpng,
                                             chartDL = charts$tmpFileDLpng,
                                             chartIS = charts$tmpFileISpng,
                                             chartIL = charts$tmpFileILpng,
                                             chartMF = charts$tmpFileMFpng,
                                             alignments = datos$alignmentdir) )
      file.rename(out, file) }
  )
  
  output$downloadMSA_clstr <- downloadHandler(
   filename = function(){
    paste(input_filename_clstrs_msa(), 'pdf', sep = '.')
   },
   content = function(file) {
    tmpFile <- tempfile(pattern = 'msa', tmpdir = '.', fileext = '.pdf')
    msaPrettyPrint(clustersMSA(), file = tmpFile, output = "pdf",
                          consensusColor="ColdHot", askForOverwrite = FALSE, 
                          showLogo="top", showLegend = FALSE,
                          subset = prettyprintCrashpreventor(),
                          furtherCode=c("\\textbf{Multiple Sequence Alignment: ClustalW} ",
                                        "\\DNAgroups{GAR,CTY}","\\defconsensus{.}{lower}{upper}",
                                        "\\showruler{1}{top}"))
    file.rename(tmpFile, file)}
   )
  
  callModule(downloadCSV, 'downloadFile2',Tabla, input_filename_generic)
  
  #CallModule only declares itself once, check all variables that should be 
  #updating automatically (reactive variables)
  callModule(downloadPAIR, 'downloadPAIR1', aln2PAIR, input_filename_pair)
  
  callModule(downloadFASTA, 'downloadFASTA1', 
             aln2FASTA, input_filename_group)
  
  callModule(downloadFASTA, 'downloadFASTAS',
             namedStringSet, input_filename_results, rows_selected)
  
  callModule(downloadFASTA, 'downloadFASTAS_clstr',
             clustersStringSet, input_filename_clstrs)
  
  session$onSessionEnded(function() {

   #Delete dir with data on session exit to avoid accumulating files
   system(paste('rm -r ', CompletePATH,'/',session$token, sep = ''))
   system('rm *.pdf')
   system('rm *.tex')
   })
