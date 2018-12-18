#---Lazy Panel Filter SERVER----


# serverMosaic <- function(input, output, session) {
  URLdata <- reactiveValues()
  DataTable <- reactiveValues()
  
  output$Unf_file <- renderDT({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    DataTable$UnfilteredDF <- read.csv(input$file1$datapath,
                                       header = input$header,
                                       sep = input$sep,
                                       quote = input$quote)
    
    if(input$disp == "head") {
      return(head(DataTable$UnfilteredDF))
    }
    else {
      return(DataTable$UnfilteredDF)
    }
    
    
  })
  
  observeEvent(input$SidebarDisp,{ if (input$SidebarDisp == 'Panel'){
    updateTabsetPanel(session ,'MainDisp', selected = 'Panel')}
    else if(input$SidebarDisp == 'Unf_file'){
      updateTabsetPanel(session ,'MainDisp', selected = 'Unf_file')
      
    }
  })
  
  output$Genes <- renderDT({
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file2)
    
    if (substring(text = input$file2$name, nchar(input$file2$name)-4,
                  nchar(input$file2$name))=='.xlsx'){
      DataTable$PanelDF <- read.xlsx(input$file2$datapath,1,
                                     header = input$header2)}
    
    if(substring(text = input$file2$name, nchar(input$file2$name)-3,
                 nchar(input$file2$name))=='.csv'){
      DataTable$PanelDF <- read.csv(input$file2$datapath,
                                    header = input$header2,
                                    sep = input$sep2,
                                    quote = input$quote2)}
    
    
    if(input$disp2 == "head") {
      return(head(DataTable$PanelDF))
    }
    else {
      return(DataTable$PanelDF)
    }
  })
  
  #filename <- reactive({paste(substring(text = input$file1$name, first = 1, 
  #last = nchar(input$file1$name)-3), 'filtered', '.csv', sep='')})
  filename <- reactive({paste(substring(text = input$file1$name, first = 1, 
                                        last = nchar(input$file1$name)-3), 
                              'filtered', sep='')})
  
  observeEvent(input$Load_filtered_file, {
    
    for (i in 1:length(unlist(DataTable$PanelDF))){
      dt_temp = filter(DataTable$UnfilteredDF, 
                       str_detect(DataTable$UnfilteredDF[,'RegionID'],
                                  as.character(DataTable$PanelDF[i,])))
      
      if (i == 1){ dt_filtered = dt_temp
      }else{
        dt_filtered = bind_rows(dt_filtered, dt_temp)
      }
    }
    DataTable$FilteredDT <- dt_filtered
    
    output$Fil_file = renderDT({
      
      if(input$disp == "head") {
        return(head(DataTable$FilteredDT))
      }
      else {
        return(DataTable$FilteredDT)
      }}
      
      ,selection = 'single')
    
  })
  #Generates hidden download button
  output$react_download_button <- renderUI({ req(input$Load_filtered_file)
    br()
    downloadButtonModule(id = 'downloadFile1', 
                         'Download CSV')})
  
  #Calls personal module for CSV file generation from DataTable object, 
  #specific inputs and filename.
  callModule(downloadCSV, 
             id = 'downloadFile1', DataTable$FilteredDT, filename)
  
  output$chooseOrg <- renderUI({
    req(input$Fil_file_rows_selected)
    wellPanel(selectInput('Organism',
                          label = 'Select Organism:',
                          choices = c('Human', 'Mus musculus')))})
  
  #Database list to choose from depends upon the organism chosen.
  output$chooseDB <- renderUI({req(input$Organism)
    wellPanel(
      if (input$Organism == 'Human'){
        selectInput('DB',
                    label = 'Select Database:',
                    choices = c('GRCH38.p12' = 'GCF_000001405.38',
                                'GRCH37.p13' = 'GCF_000001405.25',
                                'HuRef' = 'GCF_000002125.1',
                                'GMH1_1.1' = 'GCF_000306695.2'))
      }else{
        selectInput('DB',
                    label = 'Select Database:',
                    choices = c('GRCm38.p6' = 'GCF_000001635.26',
                                'Mm_Celera' = 'GCF_000002165.2'))
        
      })})
  #Reactive link generation.
  output$Link <- renderUI({req(input$Fil_file_rows_selected)
    
    wellPanel(helpText(a('Genome Data Viewer',
                         href = paste('https://www.ncbi.nlm.nih.gov/genome/gdv/?context=genome',
                                      paste('acc=', input$DB, sep = ''),
                                      URLdata$chr,
                                      URLdata$from,
                                      URLdata$to,
                                      sep='&'), 
                         target="_blank")))})
  
  output$Manual_Options <- renderUI({req(input$Fil_file_rows_selected)
    
    wellPanel(textInput('input_chr',
                        'Chromosome:',
                        value = substring(
                          DataTable$FilteredDT[input$Fil_file_rows_selected,1],
                          4,nchar(as.character(
                            DataTable$FilteredDT[input$Fil_file_rows_selected,1])
                          ))),
              
              wellPanel(numericInput('input_from',
                                     'Gap Start:',
                                     min = 1,
                                     value = DataTable$FilteredDT[input$Fil_file_rows_selected,2])),
              wellPanel(numericInput('input_to',
                                     'Gap End:',
                                     min = 1,
                                     value = DataTable$FilteredDT[input$Fil_file_rows_selected,3]))
    )})
  
  observeEvent(input$input_chr, {URLdata$chr <- paste('chr=', input$input_chr, 
                                                      sep='')})
  observeEvent(input$input_from, {URLdata$from <- paste('from=', 
                                                        input$input_from, 
                                                        sep='')})
  observeEvent(input$input_to, {URLdata$to <- paste('to=', input$input_to, 
                                                    sep='')})
# }
