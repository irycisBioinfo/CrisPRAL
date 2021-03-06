
  # set up reactives
  values <- reactiveValues()
  
  ## wait timer
  autoInvalidate <- reactiveTimer(2000)
  volumes <- c(Home = fs::path_home())

  fastqcInput <- reactive({
   
   req(input$FastqcinputFiles)
   
   # values$FastqfileSelected <- parseFilePaths(volumes,
   #                                input$FastqcinputFiles)
   withProgress(message = 'Initializing analysis', {
    values$fastqcDir <- file.dir(input$FastqcinputFiles$datapath[1])
    
    for(i in 1:length(input$FastqcinputFiles$datapath)){
     system(paste('mv ', input$FastqcinputFiles$datapath[i],' ' ,
                  values$fastqcDir,'/',input$FastqcinputFiles$name[i], sep = '' ))
    }
    
    system(paste('fastqc ', values$fastqcDir, '/*' ,sep = ''))
    
    fileList <- list.files(values$fastqcDir, pattern = "fastqc.zip", full.names = TRUE)
    fileListNames <- list.files(values$fastqcDir, pattern = "fastqc.zip", full.names = FALSE)
    
    fdl <- ngsReports::FastqcDataList(fileList)
    
    
    system(paste('cd ', values$fastqcDir, ' && ','tar -czf ', 'fastqc_reports.tgz ', 
                 str_c(fileListNames, collapse = ' '), ' && cd ', CompletePATH, sep = ''))
    
    # 
    # for(i in 1:length(input$FastqcinputFiles$datapath)){
    #  fdl@.Data[[i]]@Summary$Filename <- rep(input$FastqcinputFiles$name[i], 12)}
    
    return(fdl)
   })
  })
  
  filename <- reactive({file = paste(values$fastqcDir,'/','fastqc_reports', sep = '')})
  callModule(downloadCompressedFile, id = 'fastqcReport' ,name = filename)
  
  data <- reactive({
   
      req(input$FastqcinputFiles)
    
    data <- fastqcInput()
    #check that input$files is empty (ie no files are selected)
    if(!length(fastqcInput())){
      #choose the files from the root directory of the current volume and show only show zip files in the loadout
      shinyFileChoose(
        input,
        "files",
        roots = fastqcDir,
        session = session,
        filetypes = "zip"
      )
      #get the metadata for the file that is selected, files is the element bade in the body script
      fileSelected <- parseFilePaths(fastqcDir,
                                     input$files)
      # make the selected file a character vector
      data <- as.character(fileSelected$datapath)
      # import the selected file(s)
      #selectedData <- FastqcDataList(fileSelected)
    }
    else{
      #if a character string is provided, it will import the data for you
      if (class(fastqcInput()) == "character"){
        FastqcDataList(data)} else data
    }
    
  })
  
  
  
  #render the UI for gcTheoretical
  output$sequencedSpecies <- renderUI({
    selectInput(
      "omicSpecies",
      "Select species",
      choices = gcAvail(ngsReports::gcTheoretical, type = input$omicsType)$Name,
      selected = "Hsapiens"
    )
    
  })
  
  
  species <- observeEvent(input$omicSpecies, {
    values$omicSpecies <- input$omicSpecies
  })
  
  #export overrepresented
  
  expOS <- reactive({
    volumes <- shinyFiles::getVolumes()
    shinyFiles::shinyDirChoose(input, "dirOS",
                               roots = volumes, session = session)
    dirSelected <- shinyFiles::parseDirPath(volumes, input$dirOS)
    as.character(dirSelected)
  })
  
  observe({
    if (length(expOS())) {
      exportOverrepresented(
        data(),
        path = paste0(expOS(), "/OverrepSequences", "-", Sys.Date()),
        n = 10,
        noAdapters = TRUE
      )
    }
  })
  
  
  #### dynamic tabs to show pass warn Fail ########
  
  # render the menu item for the Summary tab
  
  output$BQflag <- renderMenu({
   
      req(input$FastqcinputFiles)
   
    if (!is.null(fastqcInput()) | length(input$files) > 1) {
      flags <- getSummary(data())
      Category <- c()
      flags <-
        subset(flags, Category == "Per base sequence quality")
      
      items <- menuItemLogic(flags = flags)
      
      values$BQflag <- items[[1]]
      values$BQcolour <- items[[2]]
      values$BQcountF <- items[[3]]
      values$BQcountW <- items[[4]]
      values$BQcountP <- items[[5]]
      
      menuItem(
        text = "Base Quality",
        tabName = "BQ",
        badgeLabel = values$BQflag,
        badgeColor = values$BQcolour
      )
      
    }
    else{
      menuItem(text = "Base Quality", tabName = "BQ")
    }
  })
  
  
  output$SQflag <- renderMenu({
   
      req(input$FastqcinputFiles)
   
    if (!is.null(fastqcInput()) | length(input$files) > 1) {
      flags <- getSummary(data())
      flags <-
        subset(flags, Category == "Per sequence quality scores")
      
      items <- menuItemLogic(flags = flags)
      
      values$SQflag <- items[[1]]
      values$SQcolour <- items[[2]]
      values$SQcountF <- items[[3]]
      values$SQcountW <- items[[4]]
      values$SQcountP <- items[[5]]
      
      menuItem(
        text = "Sequence Quality",
        tabName = "SQ",
        badgeLabel = values$SQflag,
        badgeColor = values$SQcolour
      )
      
    }
    else{
      menuItem(text = "Sequence Quality", tabName = "SQ")
    }
  })
  
  output$SCflag <- renderMenu({
   
      req(input$FastqcinputFiles)
   
    if (!is.null(fastqcInput()) | length(input$files) > 1) {
      flags <- getSummary(data())
      flags <-
        subset(flags, Category == "Per base sequence content")
      
      items <- menuItemLogic(flags = flags)
      
      values$SCflag <- items[[1]]
      values$SCcolour <- items[[2]]
      values$SCcountF <- items[[3]]
      values$SCcountW <- items[[4]]
      values$SCcountP <- items[[5]]
      
      menuItem(
        text = "Sequence Content",
        tabName = "SC",
        badgeLabel = values$SCflag,
        badgeColor = values$SCcolour
      )
      
    }
    else{
      menuItem(text = "Sequence Content", tabName = "SC")
    }
  })
  
  output$GCflag <- renderMenu({
   
      req(input$FastqcinputFiles)
   
    if (!is.null(fastqcInput()) | length(input$files) > 1) {
      flags <- getSummary(data())
      flags <-
        subset(flags, Category == "Per sequence GC content")
      
      items <- menuItemLogic(flags = flags)
      
      values$GCflag <- items[[1]]
      values$GCcolour <- items[[2]]
      values$GCcountF <- items[[3]]
      values$GCcountW <- items[[4]]
      values$GCcountP <- items[[5]]
      
      menuItem(
        text = "GC Content",
        tabName = "GC",
        badgeLabel = values$GCflag,
        badgeColor = values$GCcolour
      )
      
    }
    else{
      menuItem(text = "GC Content", tabName = "GC")
    }
  })
  
  output$NCflag <- renderMenu({
   
      req(input$FastqcinputFiles)
   
    if (!is.null(fastqcInput()) | length(input$files) > 1) {
      flags <- getSummary(data())
      flags <- subset(flags, Category == "Per base N content")
      
      items <- menuItemLogic(flags = flags)
      
      values$NCflag <- items[[1]]
      values$NCcolour <- items[[2]]
      values$NCcountF <- items[[3]]
      values$NCcountW <- items[[4]]
      values$NCcountP <- items[[5]]
      
      menuItem(
        text = "N Content",
        tabName = "NC",
        badgeLabel = values$NCflag,
        badgeColor = values$NCcolour
      )
      
    }
    else{
      menuItem(text = "N Content", tabName = "NC")
    }
  })
  
  output$SLDflag <- renderMenu({
   
      req(input$FastqcinputFiles)
   
    if (!is.null(fastqcInput()) | length(input$files) > 1) {
      flags <- getSummary(data())
      flags <-
        subset(flags, Category == "Sequence Length Distribution")
      
      items <- menuItemLogic(flags = flags)
      
      values$SLDflag <- items[[1]]
      values$SLDcolour <- items[[2]]
      values$SLDcountF <- items[[3]]
      values$SLDcountW <- items[[4]]
      values$SLDcountP <- items[[5]]
      
      menuItem(
        text = "Sequence Length Distribution",
        tabName = "SLD",
        badgeLabel = values$SLDflag,
        badgeColor = values$SLDcolour
      )
      
    }
    else{
      menuItem(text = "Sequence Length Distribution", tabName = "SLD")
    }
  })
  
  output$SDLflag <- renderMenu({
   
      req(input$FastqcinputFiles)
   
    if (!is.null(fastqcInput()) | length(input$files) > 1) {
      flags <- getSummary(data())
      flags <-
        subset(flags, Category == "Sequence Duplication Levels")
      
      items <- menuItemLogic(flags = flags)
      
      values$SDLflag <- items[[1]]
      values$SDLcolour <- items[[2]]
      values$SDLcountF <- items[[3]]
      values$SDLcountW <- items[[4]]
      values$SDLcountP <- items[[5]]
      
      menuItem(
        text = "Sequence Duplication Levels",
        tabName = "SDL",
        badgeLabel = values$SDLflag,
        badgeColor = values$SDLcolour
      )
      
    }
    else{
      menuItem(text = "Sequence Duplicaiton Levels", tabName = "SDL")
    }
  })
  
  output$OSflag <- renderMenu({
   
      req(input$FastqcinputFiles)
   
    if (!is.null(fastqcInput()) | length(input$files) > 1) {
      flags <- getSummary(data())
      flags <-
        subset(flags, Category == "Overrepresented sequences")
      
      items <- menuItemLogic(flags = flags)
      
      values$OSflag <- items[[1]]
      values$OScolour <- items[[2]]
      values$OScountF <- items[[3]]
      values$OScountW <- items[[4]]
      values$OScountP <- items[[5]]
      
      menuItem(
        text = "Overrepresented Sequences",
        tabName = "OS",
        badgeLabel = values$OSflag,
        badgeColor = values$OScolour
      )
      
    }
    else{
      menuItem(text = "Overrepresented Sequences", tabName = "OS")
    }
  })
  
  
  output$ACflag <- renderMenu({
   
      req(input$FastqcinputFiles)
   
    if (!is.null(fastqcInput()) | length(input$files) > 1) {
      flags <- getSummary(data())
      flags <- subset(flags, Category == "Adapter Content")
      
      items <- menuItemLogic(flags = flags)
      
      values$ACflag <- items[[1]]
      values$ACcolour <- items[[2]]
      values$ACcountF <- items[[3]]
      values$ACcountW <- items[[4]]
      values$ACcountP <- items[[5]]
      
      menuItem(
        text = "Adapter Content",
        tabName = "AC",
        badgeLabel = values$ACflag,
        badgeColor = values$ACcolour
      )
      
    }
    else{
      menuItem(text = "Adapter Content", tabName = "AC")
    }
  })
  
  output$KCflag <- renderMenu({
   
      req(input$FastqcinputFiles)
   
    if (!is.null(fastqcInput()) | length(input$files) > 1) {
      flags <- getSummary(data())
      flags <- subset(flags, Category == "Kmer Content")
      
      items <- menuItemLogic(flags = flags)
      
      values$KCflag <- items[[1]]
      values$KCcolour <- items[[2]]
      values$KCcountF <- items[[3]]
      values$KCcountW <- items[[4]]
      values$KCcountP <- items[[5]]
      
      menuItem(
        text = "K-mer Content",
        tabName = "KC",
        badgeLabel = values$KCflag,
        badgeColor = values$KCcolour
      )
      
    }
    else{
      menuItem(text = "K-mer Content", tabName = "KC")
    }
  })
  
  observe({
    input$files
    
    output$BQboxF <-
      renderValBox(
        count = values$BQcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
    
    output$SQboxF <-
      renderValBox(
        count = values$SQcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
    
    output$SCboxF <-
      renderValBox(
        count = values$SCcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
    
    output$GCboxF <-
      renderValBox(
        count = values$GCcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
    
    output$NCboxF <-
      renderValBox(
        count = values$NCcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
    
    output$SLDboxF <-
      renderValBox(
        count = values$SLDcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
    
    output$SDLboxF <-
      renderValBox(
        count = values$SDLcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
    
    output$OSboxF <-
      renderValBox(
        count = values$OScountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
    
    output$ACboxF <-
      renderValBox(
        count = values$ACcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
    
    output$KCboxF <-
      renderValBox(
        count = values$KCcountF,
        status = "FAIL",
        ic = "times",
        c = "red"
      )
    
    #render warn value boxes
    
    output$BQboxW <-
      renderValBox(
        count = values$BQcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
    
    output$SQboxW <-
      renderValBox(
        count = values$SQcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
    
    output$SCboxW <-
      renderValBox(
        count = values$SCcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
    
    output$GCboxW <-
      renderValBox(
        count = values$GCcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
    
    output$NCboxW <-
      renderValBox(
        count = values$NCcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
    
    output$SLDboxW <-
      renderValBox(
        count = values$SLDcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
    
    output$SDLboxW <-
      renderValBox(
        count = values$SDLcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
    
    output$OSboxW <-
      renderValBox(
        count = values$OScountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
    
    output$ACboxW <-
      renderValBox(
        count = values$ACcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
    
    output$KCboxW <-
      renderValBox(
        count = values$KCcountW,
        status = "WARN",
        ic = "exclamation",
        c = "yellow"
      )
    
    #render Pass value boxes
    
    output$BQboxP <-
      renderValBox(
        count = values$BQcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
    
    output$SQboxP <-
      renderValBox(
        count = values$SQcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
    
    output$SCboxP <-
      renderValBox(
        count = values$SCcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
    
    output$GCboxP <-
      renderValBox(
        count = values$GCcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
    
    output$NCboxP <-
      renderValBox(
        count = values$NCcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
    
    output$SLDboxP <-
      renderValBox(
        count = values$SLDcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
    
    output$SDLboxP <-
      renderValBox(
        count = values$SDLcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
    
    output$OSboxP <-
      renderValBox(
        count = values$OScountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
    
    output$ACboxP <-
      renderValBox(
        count = values$ACcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
    
    output$KCboxP <-
      renderValBox(
        count = values$KCcountP,
        status = "PASS",
        ic = "check",
        c = "green"
      )
  })
  #render fail value boxes
  
  
  
  
  
  
  # render plots
  
  ####################
  # Summary
  ####################
  
  output$SummaryFlags <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      plotSummary(
        data(),
        usePlotly = TRUE,
        cluster = input$Sumcluster,
        dendrogram = input$Sumcluster
      ) %>%
        layout(margin = list(r = 200))
    }
  })
  
  ####################
  # Read Totals
  ####################
  
  output$ReadTotals <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      plotReadTotals(data(),
                     usePlotly = TRUE,
                     duplicated = input$showDup) %>%
        layout(margin = list(l = 100, r = 200))
    }
  })
  
  ####################
  # Base Quality
  ####################
  
  output$baseQualHeatmap <- renderPlotly({
       req(input$FastqcinputFiles)
   
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      plotBaseQuals(
        data(),
        usePlotly = TRUE,
        plotType = "heatmap",
        plotValue = input$BQplotValue,
        cluster = input$BQcluster,
        dendrogram = input$BQcluster
      ) %>%
        layout(margin = list(r = 200))
    }
  })
  
  output$BaseQualitiesSingle <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      if (is.null(event_data("plotly_click")$key[[1]])) {
        num <- 1
      } else {
        click <- event_data("plotly_click")
        
        # num <- which(fqName(FastqcDataList(data())) == click$key[[1]])
        num <- which(fqName(data()) == click$key[[1]])
      }
      sub_fdl <- data()[[num]]
      plotBaseQuals(sub_fdl, usePlotly = TRUE) %>%
        layout(margin = list(r = 200, b = 50))
    }
  })
  
  ####################
  # Sequence Quality
  ####################
  
  output$seqQualHeatmap <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      plotSeqQuals(
        data(),
        cluster = input$SQcluster,
        counts = input$SQType == "Counts",
        dendrogram = input$SQcluster,
        usePlotly = TRUE
      ) %>% layout(margin = list(r = 200))
    }
  })
  
  output$SeqQualitiesSingle <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      if (is.null(event_data("plotly_click")$key[[1]])) {
        num <- 1
      } else {
        click <- event_data("plotly_click")
        # num <- which(fqName(FastqcDataList(data())) == click$key[[1]])
        num <- which(fqName(data()) == click$key[[1]])
      }
      sub_fdl <- data()[[num]]
      qualPlot <-
        plotSeqQuals(sub_fdl, usePlotly = TRUE) %>%
        layout(margin = list(r = 200),
               legend = list(orientation = 'h', title = ""))
    }
  })
  
  ####################
  # Sequence Content
  ####################
  
  output$SCHeatmap <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      plotSeqContent(
        data(),
        cluster = input$SCcluster,
        dendrogram = input$SCcluster,
        usePlotly = TRUE
      ) %>%
        layout(margin = list(r = 200))
    }
  })
  
  
  output$SCsingle <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      if (is.null(event_data("plotly_click")$key[[1]])) {
        num <- 1
      } else {
        click <- event_data("plotly_click")
        # num <- which(fqName(FastqcDataList(data())) == click$key[[1]])
        num <- which(fqName(data()) == click$key[[1]])
      }
      sub_fdl <- data()[[num]]
      plotSeqContent(sub_fdl, usePlotly = TRUE) %>%
        layout(margin = list(r = 200))
    }
  })
  
  
  ####################
  # GC Content
  ####################
  
  output$theoreticalGC <- renderUI({
    if (input$theoreticalGC) {
      radioButtons(
        inputId = "theoreticalType",
        label = "What type of data?",
        choices = c("Genome", "Transcriptome"),
        selected = "Genome"
      )
    }
  })
  
  output$GCspecies <- renderUI({
    if (!is.null(input$theoreticalGC)) {
      if (input$theoreticalGC) {
        selectInput(
          "GCspecies",
          "Select species",
          choices = gcAvail(ngsReports::gcTheoretical, type = input$theoreticalType)$Name,
          selected = "Hsapiens")
        
      }
    }
  })
  
  
  
  output$GCheatmap <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      if (is.null(input$GCspecies)) {
        GCspecies <- FALSE
      } else
        GCspecies <- input$GCspecies
      
      GCtype <- input$GCheatType == "Count"
      
      if (is.null(input$theoreticalGC)) {
        plotGcContent(
          data(),
          cluster = input$GCcluster,
          plotType = "heatmap",
          theoreticalType = input$theoreticalType,
          dendrogram = input$GCcluster,
          usePlotly = TRUE
        ) %>%
          layout(margin = list(r = 200))
      } else{
        plotGcContent(
          data(),
          cluster = input$GCcluster,
          plotType = "heatmap",
          theoreticalType = input$theoreticalType,
          theoreticalGC = input$theoreticalGC,
          dendrogram = input$GCcluster,
          species = GCspecies,
          usePlotly = TRUE
        ) %>%
          layout(margin = list(r = 200))
      }
      
    }
  })
  
  output$GCSingle <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      if (is.null(input$GCspecies)) {
        GCspecies <- FALSE
      } else
        GCspecies <- input$GCspecies
      
      if (is.null(event_data("plotly_click")$key[[1]])) {
        num <- 1
      } else {
        click <- event_data("plotly_click")
        # num <- which(fqName(FastqcDataList(data())) == click$key[[1]])
        num <- which(fqName(data()) == click$key[[1]])
      }
      sub_fdl <- data()[[num]]
      if (is.null(input$theoreticalGC)) {
        GCSingle <- plotGcContent(sub_fdl, usePlotly = TRUE,
                                  counts = FALSE)
      } else{
        GCSingle <- plotGcContent(
          sub_fdl,
          usePlotly = TRUE,
          counts = FALSE,
          theoreticalGC = input$theoreticalGC,
          theoreticalType = input$theoreticalType,
          species = GCspecies
        )
      }
      GCSingle %>%
        layout(margin = list(r = 200),
               legend = list(orientation = 'h', title = ""))
    }
  })
  
  
  ####################
  # N Content
  ####################
  
  output$NCheatmap <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      plotNContent(
        data(),
        cluster = input$Ncluster,
        dendrogram = input$Ncluster,
        usePlotly = TRUE
      )
    }
  })
  
  
  # N Content single plot
  
  output$NCsingle <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      if (is.null(event_data("plotly_click")$key[[1]])) {
        num <- 1
      } else {
        click <- event_data("plotly_click")
        #num <- which(fqName(FastqcDataList(data())) == click$key[[1]])
        num <- which(fqName(data()) == click$key[[1]])
      }
      sub_fdl <- data()[[num]]
      plotNContent(sub_fdl, usePlotly = TRUE) %>%
        layout(margin = list(r = 200, b = 50))
    }
  })
  
  
  ####################
  # Sequence Length Distribution
  ####################
  
  output$SLHeatmap <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      plotSeqLengthDistn(
        data(),
        cluster = input$SLcluster,
        dendrogram = input$SLcluster,
        counts = input$SLType == "Counts",
        usePlotly = TRUE
      ) %>%
        layout(margin = list(r = 200))
    }
  })
  
  
  output$SLSingle <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      if (is.null(event_data("plotly_click")$key[[1]])) {
        num <- 1
      } else {
        click <- event_data("plotly_click")
        # num <- which(fqName(FastqcDataList(data())) == click$key[[1]])
        num <- which(fqName(data()) == click$key[[1]])
      }
      sub_fdl <- data()[[num]]
      plotSeqLengthDistn(sub_fdl,
                         usePlotly = TRUE,
                         plotType = "line") %>%
        layout(margin = list(r = 200))
    }
  })
  
  ####################
  # Sequence Duplicaiton Levels
  ####################
  
  output$DupHeatmap <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      plotDupLevels(
        data(),
        cluster = input$Dupcluster,
        dendrogram = input$Dupcluster,
        usePlotly = TRUE
      ) %>%
        layout(margin = list(r = 200))
    }
  })
  
  output$DupSingle <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      if (is.null(event_data("plotly_click")$key[[1]])) {
        num <- 1
      } else {
        click <- event_data("plotly_click")
        # num <- which(fqName(FastqcDataList(data())) == click$key[[1]])
        num <- which(fqName(data()) == click$key[[1]])
      }
      sub_fdl <- data()[[num]]
      plotDupLevels(sub_fdl, usePlotly = TRUE) %>%
        layout(margin = list(r = 200))
    }
  })
  
  ####################
  # Overrepresented sequences
  ####################
  
  output$OSummary <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      plotOverrep(
        data(),
        usePlotly = TRUE,
        cluster = input$OScluster,
        dendrogram = input$OScluster
      ) %>%
        layout(margin = list(r = 200))
    }
    
  })
  
  
  output$OSsingle <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      if (is.null(event_data("plotly_click")$key[[1]])) {
        num <- 1
      } else {
        click <- event_data("plotly_click")
        # num <- which(fqName(FastqcDataList(data())) == click$key[[1]])
        num <- which(fqName(data()) == click$key[[1]])
      }
      sub_fdl <- data()[[num]]
      plotOverrep(sub_fdl, usePlotly = TRUE) %>%
        layout(margin = list(r = 200))
    }
  })
  
  ####################
  # Adapter Content
  ####################
  
  
  output$ACheatmap <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      ACplot <- plotAdapterContent(
        data(),
        adapterType = input$ACtype,
        usePlotly = TRUE,
        dendrogram = input$ACcluster,
        cluster = input$ACcluster
      )
      if (!is.null(ACplot))
        ACplot %>% layout(margin = list(r = 200))
      else
        stop(
          paste(
            "Sequences did not contain any",
            input$ACtype,
            "content, Please load another."
          )
        )
    }
  })
  
  output$ACsingle <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      if (is.null(event_data("plotly_click")$key[[1]])) {
        num <- 1
      } else {
        click <- event_data("plotly_click")
        # num <- which(fqName(FastqcDataList(data())) == click$key[[1]])
        num <- which(fqName(data()) == click$key[[1]])
      }
      sub_fdl <- data()[[num]]
      ACsing <- plotAdapterContent(sub_fdl, usePlotly = TRUE)
      
      if (!is.null(ACsing))
        ACsing %>%
        layout(margin = list(r = 200),
               legend = list(orientation = 'h', title = ""))
      else
        stop(
          paste(
            "Sequences did not contain any",
            input$ACtype,
            "content, Please load another."
          )
        )
    }
  })
  
  
  ####################
  # k-mer Content
  ####################
  
  output$Kheatmap <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      Kplot <- plotKmers(
        data(),
        usePlotly = TRUE,
        cluster = input$KMcluster,
        dendrogram = input$KMcluster
      )
      if (!is.null(Kplot))
        Kplot %>% layout(margin = list(r = 200))
      else
        stop(paste("Samples have no Kmer content"))
    }
  })
  
  output$Ksingle <- renderPlotly({
    if (!length(data())) {
      stop("Please load data to display plot.")
    }
    else{
      if (is.null(event_data("plotly_click")$key[[1]])) {
        num <- 1
      } else {
        click <- event_data("plotly_click")
        # num <- which(fqName(FastqcDataList(data())) == click$key[[1]])
        num <- which(grepl(click$key[[1]],
                           fqName(data())))
      }
      sub_fdl <- data()[[num]]
      Ksing <- plotKmers(sub_fdl, usePlotly = TRUE)
      
      if (!is.null(Ksing))
        Ksing %>%
        layout(margin = list(r = 200, b = 50))
      else
        stop(paste(
          "Library did not contain any identified
            Kmers please load another."
        ))
    }
  })
# }