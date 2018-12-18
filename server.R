#All For One main server - Combines mosaic analyser and Lazy panel filtering.

#First 740 lines correspond to the Mosaic Analyser.

#Lazy Panel Filtering starts at line 748.

# Big Issue: webshot does not work in Linux Debian, must comment out all graph
# printing for pdf output.

shinyServer(function(input, output, session) {
  
  #Uncomment for app to stop when session is closed:
  
  # session$onSessionEnded(function() {
  #   stopApp()
  # })
  
  #---Mosaic Analyser SERVER----
  #---Reactive values----  
  options(shiny.maxRequestSize = 100 * 1024 ^ 2)
  
  datos <- reactiveValues()
  charts <- reactiveValues()
  file_upload <- reactiveValues(upload_Adapter_state = NULL)
  file_upload <- reactiveValues(upload_Primer_state = NULL)
  
  rows_selected <- reactive({ if( length(input$tablaD_rows_selected ) != 0 ){ 
    return(input$tablaD_rows_selected )
  }else if ( input$all_clusters ){ 
    return(c(1:length(datos$Sequences)))
  }else{ return(c(1:10)) }})
  
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
      T_Refence_dir = paste(str_sub(input$Reference$datapath, 1,-8),
                            'Filtered_Reference.fasta', sep = '')
      
      system(paste(CompletePATH, '/bin/cutadapt/cutadapt -n 2 -g ', 
                   as.character(PrimerR1_5()), ' -a ', 
                   as.character(PrimerR1_3()), ' -o ', 
                   str_sub(input$Reference$datapath, 1,-8), 
                   'Filtered_Reference.fasta ', 
                   input$Reference$datapath, sep = '' ))
      return(T_Refence_dir)
    }
    else{
      return(input$Reference$datapath)
    }
  })
  
  Target_FileInput <- callModule(inputReset_server, 'reset_Target')
  Target <- reactive({
    if(file_upload$upload_Primer_state == 'uploaded' && !is.null(Target_FileInput())){ 
      T_Target_dir = paste(str_sub(Target_FileInput()$datapath, 1,-8),'Filtered_Target.fasta', sep = '')
      system(paste(CompletePATH, '/bin/cutadapt/cutadapt -n 2 -g ', as.character(PrimerR1_5()), ' -a ', as.character(PrimerR1_3()),
                   ' -o ', T_Target_dir, ' ', Target_FileInput()$datapath, sep = '' ))
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
  
  output$is.nullTarget1 <- renderUI({#Generates reactive UI depending on the presence of a Target or not.
    selectInput(
      "alnTo",
      "Alignment to:",
      if (!is.null(Target_FileInput())){c("Reference","Target")}else{c("Reference")})
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
      if(input$input_type_P == 'File_input_P'){return(readDNAStringSet(Primer_FileInput()$datapath))
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
  
  #----Processing----
  
  observeEvent(input$Accept, {
    
    req(input$R1, input$Reference)
    
    withProgress(message = "Performing analisis: ", 
                 detail = 'Clustering, may take a while...', {
      
      command = paste(
        CompletePATH,
        "/pipeline_v6.pl --r1 ",   ### Automatizar el PATH del script -- Doned
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
        collapse = "",
        sep = ""
      )
      
      system(command)#-Executes perl script
      
      incProgress( 1/5 ,detail = 'Alignments')
      
      datos$fasta = read_tsv("cluster.tsv", col_names = FALSE)
      datos$cluster = read_tsv("cluster.bak.tsv", col_names = FALSE)
      
      if (isEmpty(colnames(datos$fasta))){#-If dataset is completely filtered 
        #due to any restriction app crashes
        
        stopApp(returnValue = 'Pipe was broken')
      }else{
        colnames(datos$fasta) = c("ID", "Sequence")
        colnames(datos$cluster) = c("ClusterN", "Length", "ID", "Identity")
      } #WE NEED TO BE ABLE TO SUPPLY ERROR INFORMATION ASAP
      
      datos$Tabla_raw = datos$cluster %>% 
        group_by(ClusterN) %>% 
        mutate(Abundance = n()) %>% 
        ungroup() %>% 
        mutate(Freq = 100 *Abundance / n()) %>%
        filter(Identity == "*") %>%
        select(ID, Abundance, Freq) %>% 
        arrange(desc(Abundance))
      
      Total_Abundance = sum(datos$Tabla_raw['Abundance'])
      
      datos$Ref = readDNAStringSet(Reference())
      
      if(!is.null(Target_FileInput())){datos$Target = readDNAStringSet(Target())}
      
      datos$Sequences <- readDNAStringSet("cluster")
      
      #Alignment fo Sequences with Reference
      ls_ClusterAln_Ref = Clusters_Alignments(datos$Sequences, datos$Ref) 
      #Function in Functios.R module
      
      incProgress( 1/8 ,detail = 'Alignments... renaming')
      
      datos$aln = ls_ClusterAln_Ref$aln
      
      incProgress( 1/8 ,detail = 'Alignments')
      
      system("rm R1_* R2_* cluster* join* good* bad* final*"); #To prevent 
      #constant cashing in on the same data #26/11 removed .cluster removal for 
      #debugging
      
      datos$Tabla = inner_join(datos$Tabla_raw, ls_ClusterAln_Ref$tmp) %>% 
        mutate(score = round(score,1), Freq = signif(Freq,2)) %>% 
        arrange(desc(Abundance))
      
      output$tablaR <- renderDT(datos$Tabla, selection = "single")
      output$tablaD <- renderDT(datos$Tabla) #Table with multiple selection for 
      #PDF formatting
      output$Print <- renderPrint(datos$Tabla[sort.default(rows_selected()),])
      output$Print2 <- renderText(paste("Total amount of Reads: ", 
                                        Total_Abundance, sep = ''))
      
      #Data with named groups for FASTAS downloads.
      datos$namedSequences <- datos$Tabla %>% 
        left_join(datos$fasta) %>% select(ID, Sequence) %>%
        mutate(ID = paste(ID, '_Group', which(datos$Tabla == ID) , sep = ''))
      datos$namedStringSet <- DNAStringSet(datos$namedSequences$Sequence)
      names(datos$namedStringSet) <- datos$namedSequences$ID
      
      if(!is.null(Target_FileInput())){ 
        
        ls_ClusterAln_Tar <- Clusters_Alignments(datos$Sequences, datos$Target)
        datos$Tabla_Target <- inner_join(datos$Tabla_raw, ls_ClusterAln_Tar$tmp) %>% 
          mutate(score = round(score,1), Freq = signif(Freq,2)) %>% 
          arrange(desc(Abundance))
        datos$alnTarget <- ls_ClusterAln_Tar$aln
        Target_loc <- which.max(datos$Tabla_Target$score)
        output$Target_Location <- renderText(paste('Target is found in Cluster:', 
                                                   Target_loc, sep = ' '))
        output$tablaT <- renderDT(datos$Tabla_Target, selection = "single")
        
      }else{
        output$Target_Location <- renderText(paste('No Target was introduced', 
                                                   sep = ' '))
      }
      
      incProgress( 1/8 ,message = 'Plotting...', detail = 'Pie chart')
      
      #-Pie-Chart-processing----
      
      Other <- datos$Tabla %>% filter(Freq < 1)
      datos$Pie_data <- select(datos$Tabla, ID, Freq) %>% filter(Freq > 1) %>% 
        add_row(ID = paste('<1% groups', 
                           '(', count(Other['Freq']) , ')', 
                           sep = ''), 
                Freq = sum(Other['Freq']))
      
      #-Pie-Chart-plotting online
      
      colors <- c('rgb(211,94,96)', 'rgb(128,133,133)', 'rgb(144,103,167)', 
                  'rgb(171,104,87)', 'rgb(114,147,203)')
      
      Groups <- Group_list(datos$Pie_data) #- Function in Functions.R module
      datos$Pie_data <- add_column(datos$Pie_data, Groups)
      
      PC_online <- plot_ly(datos$Pie_data, labels = ~Groups, values = ~Freq, 
                           type = 'pie',
                           textposition = NULL,
                           hoverinfo = 'text',
                           marker = list(colors = colors, 
                                         line = list(color = '#FFFFFF', 
                                                     width = 1)), 
                           text = ~Groups, showlegend = FALSE) %>%
        layout(title = 'Main Groups',
               xaxis = list(showgrid = FALSE, zeroline = FALSE, 
                            showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, 
                            showticklabels = FALSE))
      
      output$Pie_summary <- renderText({print(PC_online)})
      
      #-Pie-Chart for plotting in PDF  
      
      PC_pdf <- plot_ly(datos$Pie_data, labels = ~ID, values = ~Freq, type = 'pie',
                        marker = list(colors = colors, 
                                      line = list(color = '#FFFFFF', width = 1)),
                        showlegend = TRUE) %>%
        layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, 
                            showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, 
                            showticklabels = FALSE))
      
      
      # charts$tmpFilePie <- tempfile(fileext = '.png')
      # export(PC_pdf, charts$tmpFilePie) #Debian-webshot error
      
      output$Pie_summary <- renderPlotly({print(PC_online)})
      
      #-----------Setting up the Data for further plotting----------------
      
      incProgress( 1/8 ,detail = 'Mutation Frequency') #User feedback
      
      #Extracting Sequences, identificators, sizes
      VCAbundance = select(datos$Tabla, Abundance)
      VCSeq = list(as.character(datos$aln@pattern))
      VCIDs = select(datos$Tabla, ID)
      SizeList = list(nchar(datos$aln@pattern))
      
      #Bind identificator with sequence and size as to not lose them when we 
      #filter by size
      Graph_Data = bind_cols(VCIDs, SizeList, VCSeq, VCAbundance)
      
      #Name the tables
      colnames(Graph_Data) = c('ID', 'Size', 'Sequence', 'Abundance')
      
      #Leave enough margin in order to correct for coverage lax
      Filtered_Graph_Data = filter(Graph_Data, Size < 250*(2-input$cov))
      Dump = filter(Graph_Data, Size > 250*(2-input$cov))
      
      TableSize = max(Filtered_Graph_Data['Size'])
      
      #----------------Processing-------------------------
      #-Mutation Frequency-This can probably be made more efficient ( was my 
      # first try at using dplyr )
      IDsSize = (count(Filtered_Graph_Data))*TableSize
      
      LetterPerPosition = data.frame(ID = rep(Filtered_Graph_Data$ID,TableSize), 
             Position = sort(rep(1:TableSize,count(Filtered_Graph_Data))),
             Chr = substring(Filtered_Graph_Data$Sequence, 
                             sort(rep(1:TableSize,count(Filtered_Graph_Data))), 
                             sort(rep(1:TableSize,count(Filtered_Graph_Data)))),
             Abundance = rep(Filtered_Graph_Data$Abundance, TableSize))
      
      LetterPerPosition_Total = LetterPerPosition %>% 
            group_by(Position, Chr, Abundance) %>% 
            count(Chr) %>% filter(Chr != '') %>% 
            transmute(Total = Abundance*n) %>% ungroup() %>%
            select(Total, Chr, Position) %>% group_by(Position, Total) %>% 
            distinct(Chr)
      
      LetterPerPosition_Norm = LetterPerPosition_Total %>% 
            ungroup() %>% group_by(Position) %>% 
            mutate_at(vars(Total),funs(./sum(Total))) %>% 
            ungroup(Position) %>% group_by(Chr, Position) %>%
            mutate_at(vars(Total), funs(sum(Total))) %>% 
            group_by(Total, Position) %>% distinct(Chr)
      
      #-Isolating deletions
      
      setProgress(detail = 'Deletions per loci')
      
      Deletions_locations = LetterPerPosition_Total %>% filter(Chr == '-') %>% 
        group_by(Position) %>% 
        mutate(Total = sum(Total)) %>% ungroup() %>% 
        select(Position, Total) %>% distinct()
      
      #-Insert locations processing
      
      setProgress(detail = 'Insertion per loci')
      
      Insert_data = as.data.frame(indel(datos$aln)@insertion)
      DFGraph_Data <- as.data.frame(Graph_Data)
      Insert_data <- Insert_data %>% 
        mutate(Count = DFGraph_Data['Abundance'][group,]) %>% 
        select(start, end, Count)
      Insert_per_loci <- as.data.frame(c(1:max(width(datos$Sequences))))
      colnames(Insert_per_loci) <- c('Position')
      Insert_per_loci <- Unravel_Positions(Insert_per_loci, Insert_data)
      
      #-----------------Plotting--------------------------
      
      setProgress(detail = 'Printing')
      
      plot.new()
      
      #-Base Frequency
      incProgress(1/30)
      
      colorscheme = c('red', 'black', 'blue', 'green', 'purple')
      colorscheme = setNames(colorscheme, c('T', 'G', 'C', 'A', '-'))
      
      MF <- plot_ly(LetterPerPosition_Norm, x = ~Position, y = ~Total, 
                    type = 'bar', name = ~Chr, color= ~Chr , 
                    colors = colorscheme) %>%
        layout(yaxis = list(title = 'Count'), barmode = 'stack')
      output$Mut_Freq <- renderPlotly({print(MF)})
      
      # charts$tmpFileMF <- tempfile(fileext = '.pdf')
      # export(MF, charts$tmpFileMF)
      
      #-Deletions per loci
      incProgress(1/30)
      
      DL <- plot_ly(Deletions_locations, x = ~Position, y = ~Total, 
                    type = 'bar',
                    name = 'Deletions') %>%
        layout(title = 'Deletions per position' , yaxis = list(title = 'Count'))
      output$Del_loc <- renderPlotly({print(DL)})
      
      # charts$tmpFileDL <- tempfile(fileext = '.pdf')
      # export(DL, charts$tmpFileDL)
      
      #-Deletion Sizes
      incProgress(1/30)
      
      Del_Data = as.data.frame(InDel_Extraction(datos$aln, VCAbundance)[1])
      
      DS = plot_ly(Del_Data, y = ~TotalD, x = ~Deletions, type = 'bar' ) %>%
        layout(title = 'Deletions per sizes' , yaxis = list(title = 'Count'),
               xaxis = list(title = 'Size'))
      output$Del_sizes <- renderPlotly({print(DS)})
      
      # charts$tmpFileDS <- tempfile(fileext = '.pdf')
      # export(DS, charts$tmpFileDS)
      
      #-Insertion Sizes
      incProgress(1/30)
      
      In_Data = as.data.frame(InDel_Extraction(datos$aln, VCAbundance)[2])
      
      IS = plot_ly(In_Data, y = ~TotalI, x = ~Insertions, type = 'bar' ) %>%
        layout(title = 'Insertions per sizes',yaxis = list(title = 'Count'),
               xaxis = list(title = 'Size'))
      output$In_sizes <- renderPlotly({print(IS)})
      
      # charts$tmpFileIS <- tempfile(fileext = '.pdf')
      # export(IS, charts$tmpFileIS)
      
      #-Insertion per Loci
      incProgress(1/30)
      
      
      IL = plot_ly(Insert_per_loci, y = ~Total_Count, x = ~Position, 
                   type = 'bar' ) %>%
        layout(title = 'Insertions per position', yaxis = list(title = 'Count'),
               xaxis = list(range = c(1, 250)))
      output$In_loc <- renderPlotly({print(IL)})
      
      # charts$tmpFileIL <- tempfile(fileext = '.pdf')
      # export(IL, charts$tmpFileIL)
      
      #- Dumped count
      
      output$Dump = renderText(paste('Filter Clusters due to size limits: ',
                                     as.character(count(Dump)),sep = ''))
      
    })
    
  })
  #tabla_rows_selected() allow alignment selection regardless of using table 
  #calculated with Reference or Target
  tabla_rows_selected <- eventReactive(c(input$tablaR_rows_selected, 
                                         input$tablaT_rows_selected),{
                                           
                                           if(input$Alignment == 'Reference'){ 
                                             return( input$tablaR_rows_selected)
                                           }else{
                                             return( input$tablaT_rows_selected)
                                           }})
  
  observeEvent(c(input$tablaR_rows_selected, input$tablaT_rows_selected, 
                 input$alnType, input$alnTo), {
                   
                   withProgress(message = "Performing analisis", {
                     
                     req(datos$Tabla)
                     #To prevent crashing WITHOUT any error being handed to you.
                     req(tabla_rows_selected()) 
                     
                     Seq <- datos$namedStringSet[tabla_rows_selected()]
                     
                     if(input$alnTo == "Reference")
                     {
                       datos$aln2 = pairwiseAlignment(datos$Ref, Seq, 
                                                      type = input$alnType)
                     }else{
                       datos$aln2 = pairwiseAlignment(datos$Target, Seq, 
                                                      type = input$alnType)
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
                     br()
                     output$downloadPAIR_file <- renderUI({
                       br()
                       downloadButtonModule('downloadPAIR1',
                                            'Download Pairwise alignment')
                       })
                     br()
                     output$Unaligned_Fasta <- renderUI({
                       wellPanel(
                         h5('Unaligned Sequence'),
                         verbatimTextOutput('unalignedID'),
                         verbatimTextOutput('unalignedFASTA'))
                      })
                     
                     output$unalignedID <- renderText(as.character(
                       datos$namedSequences[tabla_rows_selected(),]$ID))
                     
                     output$unalignedFASTA <- renderText(as.character(
                       datos$namedSequences[tabla_rows_selected(),]$Sequence))
                   })
                 })
  
  #-Reactive functions for CallModule
  
  aln2PAIR <- reactive({datos$aln2})
  aln2FASTA <- reactive({alignedSubject(datos$aln2)})
  namedStringSet <- reactive({datos$namedStringSet})
  
  input_filename_generic <- reactive({paste(
                                      str_replace(#Removes posible double 
                                                  #underscore "_"
                                      str_replace(#Removes Read number ID
                                        substring(text = input$R1$name, 
                                                  first = 1, 
                                                  last = nchar(input$R1$name)-6)
                                        ,'R1',''), '__',''), sep='')})
  input_filename_results <- reactive({paste( input_filename_generic(), 
                                             '_Results', sep='' )})
  input_filename_group <- reactive({paste( input_filename_generic(), 
                                           '_Group', tabla_rows_selected(), 
                                           sep='')})
  input_filename_pair <- reactive({paste( input_filename_group() , 
                                          '_pairwise_alignment', sep='')})
  
  
  datos$alignmentdir <- 'empty' #Temporal asignment
  
  output$downloadReport <- downloadHandler(
    filename = function() {
      paste(input_filename_results(),'pdf', sep = '.')
    },
    content = function(file) {
      src <- normalizePath('reports.Rmd')
      out <- rmarkdown::render( 'reports.Rmd', params = list(selection = rows_selected(),
                                                             table = datos$Tabla,
                                                             chartPie = charts$tmpFilePie,
                                                             chartDS = charts$tmpFileDS,
                                                             chartDL = charts$tmpFileDL,
                                                             chartIS = charts$tmpFileIS,
                                                             chartIL = charts$tmpFileIL,
                                                             chartMF = charts$tmpFileMF,
                                                             alignments = datos$alignmentdir) )
      file.rename(out, file) }
  )
  
  callModule(downloadCSV, 'downloadFile2', 
             datos$Tabla, input_filename_generic)
  
  #CallModule only declares itself once, check all variables that should be 
  #updating automatically (reactive variables)
  
  callModule(downloadPAIR, 'downloadPAIR1', aln2PAIR, input_filename_pair)
  
  callModule(downloadFASTA, 'downloadFASTA1', 
             aln2FASTA, input_filename_group)
  
  callModule(downloadFASTA, 'downloadFASTAS',
             namedStringSet, input_filename_results, rows_selected)
  
  
  #---Lazy Panel Filter SERVER----
  
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
  
})
