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

#-Compartmentalization line-
dir.create(session$token)

# Selected entries for download will be determined by the tab selected.
download_table_selected <- reactive({
  if(input$Alignment2 == 'Reference'){
    return(input$tablaD_rows_selected)
  }else{
    return(input$tablaTD_rows_selected)
}})

rows_selected <- reactive({
  if( length( download_table_selected() ) != 0 ){ 
    return( download_table_selected() )
  }else if ( input$all_clusters ){ 
    return(c(1:length(datos$Sequences)))
  }else{ return(c(1:10)) }
})

Table_to_download <- reactive({
  if(input$Alignment2 == 'Reference'){
    return(Tabla())
  }else{
    return(TablaT())
}})

warningCLSTR <- reactive({if(length(datos$clustersStringSet) > 100){
    return('Warning: Number of sequences for Multiple Sequence Alignment exceed max computing limit, 100 will be computed but only the first 50 entries will be shown')
  }else if(length(datos$clustersStringSet) > 50){
    return('Warning: Number of sequences for Multiple Sequence Alignment exceed max printing limit, only the first 50 entries will be shown')
  }else{return(FALSE)}
})

#---Reactive values-Data storage----  

datos <- reactiveValues()
charts <- reactiveValues()
checks <- reactiveValues()
checks$Table_has_been_edited <- FALSE
checks$Table_back_to_default <- FALSE
file_upload <- reactiveValues(upload_Adapter_state = NULL)
file_upload <- reactiveValues(upload_Primer_state = NULL)

#-Reactive input to allow for compressed files processing.
R1_file <- reactive({
 
 if(checks$example){return('../basic_Trials/Mutation_mus_R1_001.fastq')}
 
if(!is.null(input$R1)){ 
 if (str_sub(input$R1$datapath, -3) == '.gz'){
  #Function "Unzip_file" defined in Functions.R
  Unzip_file(CompletePATH, input$R1$datapath)}
 else{
  input$R1$datapath}}
})
R2_file <- reactive({
 if(checks$example){return('../basic_Trials/Mutation_mus_R2_001.fastq')}
 
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
 
 if(checks$example){return('../basic_Trials/Reference.fasta')}
 
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

Target_location <- reactive({find_target_location(Tabla(), TablaT(), datos$Target)}) #function in Functions.R

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

#----Alignment and Table filtering dinamic UI----

Back_to_default <- reactive({
 if(!checks$Table_back_to_default){return(actionButton('Default', 'Default'))}
 else{return()}
})

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
           c("global", "local", "overlap", "global-local", "local-global"), 
           selected = input$general_allType,
          ),
          uiOutput('is.nullTarget1'),
          verbatimTextOutput('Print2'),
          uiOutput('coloured_text'),
          verbatimTextOutput('Target_Location'))),
  column(5, 
         wellPanel(h3("Cluster filtering"),
                   sliderInput('score_threshold',label = p('Scoring threshold'), 
                               min = floor(min(datos$Tabla_Original %>% select(score))),
                               max = round(max(datos$Tabla %>% select(score))), 
                               value = floor(min(datos$Tabla %>% select(score)))),
                   br(),
                   sliderInput('abundance_minimum',label = p('Abundance minimun'), 
                               min = floor(min(datos$Tabla_Original %>% select(Abundance))),
                               max = round(max(datos$Tabla %>% select(Abundance))), 
                               value = floor(min(datos$Tabla %>% select(Abundance)))),
                   fluidRow(
                    column(2,
                           actionButton('Filter', 'Accept')),
                    column(2, Back_to_default())
                   )
         )))})

observeEvent(c(input$tablaR_cell_edit, input$tablaT_cell_edit),{
 
 edit <- list(tablaR = input$tablaR_cell_edit, tablaT = input$tablaT_cell_edit)
 
 info = edit[edit != 'NULL']
 
 # str(info) Check info for debugging purposes
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

input_filename_csv <- reactive({paste( input_filename_generic(), 
                                       '_',input$Alignment2, sep='' )})

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

callModule(downloadCSV, 'downloadFile2',Table_to_download, input_filename_csv)

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
 # system('rm *.pdf')
 # system('rm *.tex')
})
