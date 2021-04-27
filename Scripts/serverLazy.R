#---Lazy Panel Filter SERVER----


observeEvent(input$go_to_lazyP, {
  
  updateTabsetPanel(session, 'app', selected = 'Lazy_Panel')
  
})

# serverMosaic <- function(input, output, session) {
  URLdata <- reactiveValues()
  DataTable <- reactiveValues()
  values <- reactiveValues()
  Charts <- reactiveValues()
  values$show <- "true"
  
  
  DataTable$EmptyDF <- read.csv("./Lazy_Filter_Project/Empty.csv",
                                 header=TRUE,
                                 sep = ',')
  
  #Generate text help box
  output$Select_column <- renderText({'Select column with genes identificators (Use "clean identifiers" if the column is contaminated), 
    gap start and gap stop coordinates'})
  
  output$Unf_file <- renderDT(extensions = c('Buttons','FixedColumns','Scroller','ColReorder'), 
                              options = list(dom = 'Bfrtip', 
                                             scrollX = TRUE,
                                             deferRender = TRUE, scrollY = 400, scroller = TRUE,
                                             search = list(smart = FALSE),
                                             buttons = list('copy', list(extend = 'collection',
                                                                         buttons = c('csv', 'excel'),
                                                                         text = 'Download'),
                                                            'colvis')), server = TRUE,{
    
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    
    readData <- function(file = input$file1$datapath, header = input$header, sep = input$sep, quote = input$quote, strip.white = TRUE){
      
      out <- tryCatch({
        
        read.table(file,
                   header = header,
                   sep = sep,
                   quote = quote,
                   strip.white = strip.white)
        
      },
      error = function(cond){
        
        if(suppressWarnings(str_detect(cond, "did not have"))){
          
          showModal(modalDialog(
            fluidPage(
              h1(strong("::ERROR::"),align="center"),
              hr(),
              fluidRow(column(4,offset = 4,
                              div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                                  img(src="caution-icon.png", height='150', width='150', align="middle")))),
              h2(paste("It seems the separation value selected is not appropriate, try selecting another.",sep = ''), align = 'center'),
              h3(paste("If problem persists contact system administrator,",datos$system_administrator,sep = ''), align = 'center'),
              
              easyClose = TRUE,
              footer = NULL
            )))
          
          return(read_lines(file))
          }else{return("Unknown error, contact administration at sergio.fern1994@gmail.com")}
        
      })
      return(out)
    }
    
    DataTable$UnfilteredDF <- readData(input$file1$datapath)
    
    if(input$disp == "head") {
      return(head(as.tibble(DataTable$UnfilteredDF)))
    }
    else {
      return(as.tibble(DataTable$UnfilteredDF))
    }
    
  }, rownames=FALSE)
  
  

  output$Edit_File <- renderDT(filter = 'top',
                               extensions = c('Buttons','FixedColumns','Scroller','ColReorder'), 
                               options = list(dom = 'Bfrtip', 
                                              scrollX = TRUE,
                                              deferRender = TRUE, scrollY = 400, scroller = TRUE,
                                              search = list(smart = FALSE),
                                              buttons = list('copy', list(extend = 'collection',
                                                                          buttons = c('csv', 'excel'),
                                                                          text = 'Download'),
                                                             'colvis')), server = TRUE,rowname = FALSE,{

    if (length(input$Fil_file2_columns_selected) == 1){ 
      
      DataTable$EditedDF <- DataTable$FilteredDT %>% 
        select(Identifier = input$Fil_file2_columns_selected)

    }
      else if (length(input$Fil_file2_columns_selected) == 2) {
       vars = c( Identifier = input$Fil_file2_columns_selected[1], Gap_Start = input$Fil_file2_columns_selected[2])
      DataTable$EditedDF <- DataTable$FilteredDT %>% select(vars)


    } else if (length(input$Fil_file2_columns_selected) == 3) {
      vars = c( Identifier = input$Fil_file2_columns_selected[1], 
                Gap_Start = input$Fil_file2_columns_selected[2], 
                Gap_Stop = input$Fil_file2_columns_selected[3])
      DataTable$EditedDF <- DataTable$FilteredDT %>%
        select(vars)
      
      DataTable$EditedDF_final = DataTable$EditedDF
      
      #---

      #---

    }else{return(DataTable$EmptyDF)}

    return(DataTable$EditedDF)
    
  })
  
  observeEvent(input$Clean,{
   
   req(input$Fil_file2_columns_selected)
    
    values$identifiers <- select(DataTable$FilteredDT, 'Identifier' = input$Fil_file2_columns_selected[1]) %>% pull("Identifier" )%>%
      str_split_fixed(".chr.",n=2)
    DataTable$FilteredDT = DataTable$FilteredDT %>% select(-input$Fil_file2_columns_selected[1]) %>% 
      add_column(Identifier = values$identifiers[,1])
    
  })
  
  output$Fil_file2 = renderDT(filter = 'top',
                              extensions = c('Buttons','FixedColumns','Scroller','ColReorder'), 
                              options = list(dom = 'Bfrtip', 
                                             scrollX = TRUE,
                                             deferRender = TRUE, scrollY = 400, scroller = TRUE,
                                             search = list(smart = FALSE),
                                             buttons = list('copy', list(extend = 'collection',
                                                                         buttons = c('csv', 'excel'),
                                                                         text = 'Download'),
                                                            'colvis')), server = TRUE,{
    
    if(input$disp == "head") {
      return(head(DataTable$FilteredDT))
    }
    else {
      return(DataTable$FilteredDT)
    }}
    ,selection = list(target = 'column'))
  
  
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
                                    quote = input$quote2, 
                                    strip.white = TRUE)}
    
    
    if(input$disp2 == "head") {
      return(head(DataTable$PanelDF))
    }
    else {
      return(DataTable$PanelDF)
    }
    
    datatable(selection = 'none')
    
  })
  
  filenameLazy <- reactive({paste(str_split(input$file1$name, pattern = "\\.")[[1]][1], 
                              'filtered', sep='_')})
  
  observeEvent(input$Load_filtered_file, {
    withProgress(message = 'Filtering data',{
    if (!is.null(DataTable$PanelDF)){
      for (i in 1:length(unlist(DataTable$PanelDF))){
        incProgress(amount=1/length(unlist(DataTable$PanelDF)))
        dt_temp = filter_all(DataTable$UnfilteredDF, 
                             any_vars(str_detect(., paste('(?<![:alnum:])',DataTable$PanelDF[i,], 
                                                          '(?![:alnum:])', sep = ''))))
        
        if (i == 1){ dt_filtered = dt_temp
        }else{
          dt_filtered = bind_rows(dt_filtered, dt_temp)
        }
      }
      DataTable$FilteredDT <- dt_filtered}
    
    else{
      DataTable$FilteredDT <- DataTable$UnfilteredDF
      
    }
    
    output$Fil_file = renderDT(filter = 'top',
                               extensions = c('Buttons','FixedColumns','Scroller','ColReorder'), 
                               options = list(dom = 'Bfrtip', 
                                              scrollX = TRUE,
                                              deferRender = TRUE, scrollY = 400, scroller = TRUE,
                                              colReorder = TRUE,
                                              search = list(smart = FALSE),
                                              buttons = list('copy', list(extend = 'collection',
                                                                          buttons = c('csv', 'excel'),
                                                                          text = 'Download'),
                                                             'colvis')), server = TRUE,{
      
      if(input$disp == "head") {
        return(head(DataTable$FilteredDT))
      }
      else {
        return(DataTable$FilteredDT)
      }}
      
      ,selection = 'single')
    })# Progress bar end brackets
  })
  ##################################################.
  ##### Deprecated ####.
  ##################################################.
  #Generates hidden download button
  # output$react_download_button <- renderUI({ req(input$Load_filtered_file)
  #   br()
  #   downloadButtonModule(id = 'downloadFile1', 
  #                        'Download CSV')})
  
  #Calls personal module for CSV file generation from DataTable object, 
  #specific inputs and filename.
  
  # FilteredDT <- reactive({DataTable$FilteredDT})
  # callModule(downloadCSV, 
  #            id = 'downloadFile1', FilteredDT, filenameLazy)
  
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
  observeEvent(input$Graph ,{
    
    req(DataTable$EditedDF_final)
    
    DataTable$PrintingDF <- DataTable$EditedDF %>% group_by(Identifier) %>%
      transmute(Gap_Size = Gap_Stop - Gap_Start) %>%
      mutate_at(vars(Gap_Size), funs(sum(Gap_Size)))

    DataTable$PrintingDF_reps <- DataTable$PrintingDF %>% group_by(Identifier) %>%
      summarise(Times_repeated=n()) %>% distinct()

    DataTable$PrintingDF <- right_join(DataTable$PrintingDF_reps, DataTable$PrintingDF,  by = "Identifier") %>% distinct()

    Charts$DF_chart <- plot_ly(DataTable$PrintingDF, x = ~Gap_Size, y = ~Times_repeated, 
                               text = DataTable$PrintingDF$Identifier, hoverinfo = 'text') %>%
      layout(title = 'Size vs times repeated',
             yaxis = list(zeroline = FALSE),
             xaxis = list(zeroline = FALSE))

    output$Size_vs_Depth <- renderPlotly({
      print(Charts$DF_chart)})
      }
  )

  observeEvent(input$input_chr, {URLdata$chr <- paste('chr=', input$input_chr, 
                                                      sep='')})
  observeEvent(input$input_from, {URLdata$from <- paste('from=', 
                                                        input$input_from, 
                                                        sep='')})
  observeEvent(input$input_to, {URLdata$to <- paste('to=', input$input_to, 
                                                    sep='')})