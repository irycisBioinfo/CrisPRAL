# serverMosaic_batch

observeEvent(input$go_to_mosaic_batch, {
  
  updateTabsetPanel(session, 'app', selected = 'Mosaic_Finder_batch')
  
})

 # Selects pipeline given the option "reverse_complement"
pipeline <- reactive({if(input$reverse_complement){
  return("/pipeline_vUMI.pl --r1 ")
}else{return("/pipeline_v7.pl --r1 ")}
})

# ############################################################################.
# ADAPTERS ####
##############################################################################.

#This following set of observeEvents allows the user to change his mind on 
#uploaded Adapter files.
observeEvent(c(input$AdapterFile_batch, input$uploadAdapters_batch), {
  file_upload$upload_Adapter_state_batch <- 'uploaded'
})
observeEvent(c(input$resetA_batch, input$input_type_A_batch), {
  file_upload$upload_Adapter_state_batch <- 'reset'
})

Adapter_Directinput_batch <- eventReactive(input$uploadAdapters_batch,
                                     { c(input$R1_Adapter_batch,input$R2_Adapter_batch) })
Adapter_FileInput_batch <- reactive({
  if (is.null(file_upload$upload_Adapter_state_batch)) {
    return(NULL)
  } else if (file_upload$upload_Adapter_state_batch == 'uploaded') {
    
    if (input$input_type_A_batch == 'File_input_A_batch'){return(input$AdapterFile_batch)}
    else{return(Adapter_Directinput_batch())}
  }
  else if (file_upload$upload_Adapter_state_batch == 'reset') {
    return(NULL)
  }
})

Adapter_sequences_batch <- reactive({
  
  if (is.null(Adapter_FileInput_batch())){
    return("Empty")
  }else if (!is.null(Adapter_FileInput_batch())){
    if(input$input_type_A_batch == 'File_input_A_batch'){
      return(readDNAStringSet(Adapter_FileInput_batch()$datapath))
    }else{
      return(DNAStringSet(Adapter_FileInput_batch()))}
  }})

AdapterR1_batch <- reactive({
  if (is.null(Adapter_FileInput_batch())){
    return("Empty")
  }
  return(as.character(Adapter_sequences_batch()[1]))
})
AdapterR2_batch <- reactive({
  if (is.null(Adapter_FileInput_batch())){
    return("Empty")
  }
  return(as.character(Adapter_sequences_batch()[2]))
})

command.adapter <- reactive({
  if( is.null(Adapter_FileInput_batch())){
    return(NULL)
  }else{return(str_c('-adapters ', Adapter_FileInput_batch()$datapath))}
})

output$checkR1A_batch <- renderText(paste("R1 Adapter: ",AdapterR1_batch(), sep = ''))
output$checkR2A_batch <- renderText(paste("R2 Adapter: ",AdapterR2_batch(), sep = ''))


# ############################################################################.
# PRIMERS ####
##############################################################################.

observeEvent(c(input$PrimerFile_batch, input$uploadPrimers_batch), {
  file_upload$upload_Primer_state_batch <- 'uploaded'
}) 
observeEvent(c(input$resetP_batch, input$input_type_P_batch), {
  file_upload$upload_Primer_state_batch <- 'reset'
})

Primer_Directinput_batch <- eventReactive(input$uploadPrimers_batch,{ 
  c(input$F_Primer1_batch,input$R_Primer2_batch) 
})
Primer_FileInput_batch <- reactive({
  if (is.null(file_upload$upload_Primer_state_batch)) {
    return(NULL)
  } else if (file_upload$upload_Primer_state_batch == 'uploaded') {
    
    if (input$input_type_P_batch == 'File_input_P_batch'){return(input$PrimerFile_batch)}
    else{return(Primer_Directinput_batch())}
  }
  else if (file_upload$upload_Primer_state_batch == 'reset') {
    return(NULL)
  }
})

#Primer upload through file:
#Primers loaded as reactive functions to prevent residual data lingering.
Primer5_sequences_batch <- reactive({
  
  if (is.null(Primer_FileInput_batch())){
    return("Empty")
  }else if (!is.null(Primer_FileInput_batch())){
    if(input$input_type_P_batch == 'File_input_P_batch'){return(readDNAStringSet(
      Primer_FileInput_batch()$datapath))
    }else{
      return(DNAStringSet(Primer_FileInput_batch()))}
  }})
Primer3_sequences_batch <- reactive({
  
  if (is.null(Primer_FileInput_batch())){
    return("Empty")
  }else if (!is.null(Primer_FileInput_batch())){
    return(reverseComplement(Primer5_sequences_batch()))
  }})

PrimerR1_5_batch <- reactive({
  if (is.null(Primer_FileInput_batch())){
    return("Empty")
  }
  return(as.character(Primer5_sequences_batch()[1]))
})
PrimerR1_3_batch <- reactive({
  if (is.null(Primer_FileInput_batch())){
    return("Empty")
  }
  return(as.character(Primer3_sequences_batch()[2]))
})
PrimerR2_5_batch <- reactive({
  if (is.null(Primer_FileInput_batch())){
    return("Empty")
  }
  return(as.character(Primer5_sequences_batch()[2]))
})
PrimerR2_3_batch <- reactive({
  if (is.null(Primer_FileInput_batch())){
    return("Empty")
  } 
  return(as.character(Primer3_sequences_batch()[1]))
})

command.primers <- reactive({
  if( is.null(Primer_FileInput_batch())){
    return(NULL)
  }else{return(str_c('-primers ', Primer_FileInput_batch()$datapath))}
})

output$Explicative_TextP_batch <- renderText(paste(
  'The file should contain four lines, two', 
  'lines with ">" indicating Primer ID and', 
  'another two lines with the sequences', 'themselves. Ex:',
  '>Primer1', 'ACGT', '>Primer2', 'GTAC',sep = '\n'))

output$Explicative_TextA_batch <- renderText(paste(
  'The file should contain four lines, two', 
  'lines with ">" indicating Adapter ID and', 
  'another two lines with the sequences', 'themselves. Ex:',
  '>Adapter1', 'ACGT', '>Adapter2', 'GTAC',sep = '\n'))

output$checkFP1_batch <- renderText(paste("Read 1, 5' end: ",PrimerR1_5_batch(), 
                                    sep = ''))
output$checkFP2_batch <- renderText(paste("Read 1, 3' end: ",PrimerR1_3_batch(), 
                                    sep = ''))
output$checkRP1_batch <- renderText(paste("Read 2, 5' end: ",PrimerR2_5_batch(), 
                                    sep = ''))
output$checkRP2_batch <- renderText(paste("Read 2, 3' end: ",PrimerR2_3_batch(), 
                                    sep = ''))

# ############################################################################.
# REFERENCE & TARGETS ####
##############################################################################.

Reference_batch <- reactive({
  
  if(file_upload$upload_Primer_state_batch == 'uploaded'){
    extensionR <- str_split(input$Reference_batch$datapath, '\\.')[[1]][length(str_split(input$Reference_batch$datapath, '\\.')[[1]])]
    lastIndexR <- nchar(extensionR)
    T_Refence_dir = paste0(str_sub(input$Reference_batch$datapath, 1,-(lastIndexR+3)),
                           'Filtered_Reference.fasta')
    
    if(extensionR == 'txt'){
      RefAsFasta <- paste0(str_sub(input$Reference_batch$datapath, 1,-(lastIndexR+3)) ,'RefAsFasta.fasta')
      system(paste( 'cp', input$Reference_batch$datapath, RefAsFasta, sep = ' '))
    }else if(extensionR == 'fasta'){
      RefAsFasta <- input$Reference_batch$datapath
    }else{
      return(':ERROR: Unknown extension, please make sure that the Reference file has extension .txt or .fasta')
    }
    
    system(paste0('cutadapt -n 2 -g ', 
                  as.character(PrimerR1_5_batch()), ' -a ', 
                  as.character(PrimerR1_3_batch()), ' -o ', 
                  T_Refence_dir,' ',
                  RefAsFasta))
    
    return(T_Refence_dir)
  }
  else{
    return(input$Reference_batch$datapath)
  }
})

# ##########################################################################.
# MISMATCHES FLAG ####
# ##########################################################################.

mismatches <- reactive({
  if(input$print_mismatches){
    return('-mismatches')
    }else{return(NULL)}
  })

# ##########################################################################.
# APPEND SEQUENCES FLAG ####
# ##########################################################################.

sequences <- reactive({
  if(input$print_sequences){
    return('-append_sequences')
  }else{return(NULL)}
})

# ##########################################################################.
# Scoring Threshold & Interest Range ####
# ##########################################################################.

scoring_threshold <- reactive({ # scoring_threshold = 'detect' if none is selected 
                                # and return the value provided if user specifies it
  if(input$detect_theshold){
    return('detect')
  }else{
    return(input$scoring_threshold)
    } 
})

interest_range <- reactive({
  if(input$interest_range_select){
    return(paste0(input$min_interest_range,'-',input$max_interest_range))
  }else{
    return(NULL)
  } 
})

command.interest_range <- reactive({
  if( is.null(interest_range())){
    return(NULL)
  }else{return( str_c('-indel_interest_range',
                      interest_range(), sep = ' ') )}
})

output$advanced_for_the_time_being <- renderText(paste0('For the time being advanced options are all fixed to default values of identity and coverage of 1 and a minimun length of 100 bps'))

# ##########################################################################.
# EXECUTING PIPELINE ####
# ##########################################################################.

observeEvent(input$Accept_batch,{
  
  # ##########################################################################.
  # MODAL DIALOG ASSERTION ####
  # ##########################################################################.
  
  required_values <- list("Reference sequence" = input$Reference_batch$name, 'file (zip, tar.gz) containing the data' = input$dir$name)
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
  #############################################################################.
  #### COMMAND ####
  #############################################################################.
  
  withProgress(message = "Performing analisis: ", detail = 'Executing pipeline', {
  
  datos$tmppipelinedir <- str_c(session$token, '/data', str_sub(tempfile(), start = 21))
  dir.create(datos$tmppipelinedir)
  
  # req(input$Reference_batch, input$dir)
  unziped_files <- conditional_input_decision_tree(input$dir)
  
  # if(unzipped_files == "Unrecognised input"){showModal error; setwd(CompletePath)}
  
  R1_files <- grep(unziped_files, pattern = '_R1_', value = TRUE)
  command <- str_c(ScriptsPATH, '/Mosaic_standalone_pipe.R',
                paste(' ',paste(R1_files,collapse = ' '),
                '-reference', 
                Reference_batch(),
                command.adapter(),
                command.primers(),
                '-alignment',
                input$alnType_batch,
                '-output',
                str_c(CompletePATH,'/',as.character(datos$tmppipelinedir)),
                mismatches(),
                '-score_threshold',
                scoring_threshold(),
                sequences(),
                command.interest_range(), sep = ' '), sep = '')

  incProgress( 1/2 ,detail = 'may take a while...')
  system(command)
  setwd(datos$tmppipelinedir)
  final_files <- dir('.', pattern = '*.xlsx', full.names = TRUE)
  path_to_zip_file <- str_c(CompletePATH,'/',as.character(datos$tmppipelinedir),'/','results.zip')
  incProgress( 1/4 ,detail = 'zipping files...')
  zip(path_to_zip_file,files = final_files)
  setwd(CompletePATH)
  
  ###################################.
  ##### FILE READY MODAL DIALOG ####
  ###################################.
  
  showModal(modalDialog(
    fluidPage(
      h1(strong("File ready"),align="center"),
      hr(),
      fluidRow(column(4,offset = 4,
                      div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                          img(src='ready_for_download.png', height='150', width='150', align="middle")))),
      h2(paste("Click to download file." ,sep = ''), align = 'center'),
      fluidRow(column(4,offset = 4.5,
                      callModule(downloadZIP, 'downloadZIP_batch', data = path_to_zip_file, 'results'))),
      easyClose = TRUE,
      footer = NULL
    )))
 
  })# withProgress closing bracket
 })# observerEvent end bracket

# input_filename_batch <- reactive({
#   'results'})
