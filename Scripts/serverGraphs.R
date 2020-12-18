#Graph generating script.

Check_if_graphs_are_altered <- eventReactive(input$Visualization,{
 
 if(!is.null(datos$Tabla_prev)){
  if(nrow(datos$Tabla_prev) != nrow(Tabla()) && input$Visualization == 'Graphics'){
   
   showModal(modalDialog(
    fluidPage(
     h1(strong("Warning"),align="center"),
     hr(),
     fluidRow(column(4, offset = 4,
                     div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                         img(src="caution-icon.png", height='150', width='150', align="middle")))),
     h2(paste("These graphics may not represent anymore the data in Clusterization and Alignment.
              Reload to ensure both results match.",
              sep = ''), align = 'center'),
     easyClose = TRUE,
     footer = NULL
    )))
   
  }
  return(FALSE)
 }})

observeEvent(input$Visualization,{
 Check_if_graphs_are_altered()
})

observeEvent(input$GG,{
 
 withProgress(message = "Plotting graphs: ", detail = 'Calculating', {
  
  if(!is.null(datos$Tabla_prev) && file.exists(paste(session$token,'/default_charts.Rdata', sep = '')) && 
     file.exists(session$token,'/defaults_charts_doc.Rdata', sep = '')){
   
   if(nrow(datos$Tabla_prev) != nrow(Tabla()) && nrow(Tabla()) == nrow(datos$Tabla_Original)){
    #This double conditioning while seeming redundant allows the re-processing
    # of graphs in the case saved objects are not loading correctly.
    
    load(paste(session$token,'/default_charts.Rdata', sep = ''))
    load(paste(session$token,'/defaults_charts_doc.Rdata', sep = ''))
    
    output$Mut_Freq <- renderPlotly({print(MF)})
    output$In_sizes <- renderPlotly({print(IS)})
    output$Del_sizes <- renderPlotly({print(DS)})
    output$In_loc <- renderPlotly({print(IL)})
    output$Del_loc <- renderPlotly({print(DL)})
    datos$Pie_data <- Pie_data
    output$Pie_summary <- renderPlotly({Edit_Pie_chart()})
    
    charts$tmpFilePiepng <- tmpFilePiepng
    charts$tmpFileDSpng <- tmpFileDSpng
    charts$tmpFileDLpng <- tmpFileDLpng
    charts$tmpFileISpng <- tmpFileISpng
    charts$tmpFileILpng <- tmpFileILpng
    charts$tmpFileMFpng <- tmpFileMFpng
    
    return()
    
   }}
  
  #-Pie-Chart-processing----
  
  Other <- datos$Tabla %>% filter(Freq < 1)
  datos$Pie_data <- select(Tabla(), ID, Freq) %>% filter(Freq > 1) %>% 
   add_row(ID = paste('<1% groups', 
                      '(', count(Other['Freq']) , ')', 
                      sep = ''), 
           Freq = sum(Other['Freq']))
  
  #-Pie-Chart-plotting online
  
  Groups <- Group_list(datos$Pie_data) #- Function in Functions.R module
  datos$Pie_data <- add_column(datos$Pie_data, Groups)
  Pie_data <- datos$Pie_data #required for default workspace saving
  
  output$Pie_summary <- renderText({print(Edit_Pie_chart())})
  
  #-Pie-Chart for plotting in PDF 
  
  colors_pie <- c('rgb(211,94,96)', 'rgb(128,133,133)', 'rgb(144,103,167)', 
                  'rgb(171,104,87)', 'rgb(114,147,203)')
  
  PC_pdf <- plot_ly(datos$Pie_data, labels = ~ID, values = ~Freq, type = 'pie',
                    marker = list(colors = colors_pie, 
                                  line = list(color = '#FFFFFF', width = 1)),
                    showlegend = TRUE) %>%
   layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, 
                       showticklabels = FALSE),
          yaxis = list(showgrid = FALSE, zeroline = FALSE, 
                       showticklabels = FALSE))
  
  
  # charts$tmpFilePie <- tempfile(fileext = '.png')
  
  charts$tmpFilePiehtml <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.html'), start = 16), sep = '')
  charts$tmpFilePiepng <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.png'), start = 16), sep = '')
  tmpFilePiepng <- charts$tmpFilePiepng
  htmlwidgets::saveWidget(PC_pdf, paste(CompletePATH,'/', charts$tmpFilePiehtml, sep = ''))
  webshot2::webshot(charts$tmpFilePiehtml, charts$tmpFilePiepng)
  
  output$Pie_summary <- renderPlotly({Edit_Pie_chart()})
  
  #-----------Setting up the Data for further plotting----------------
  
  incProgress( 1/8 ,detail = 'Mutation Frequency') #User feedback
  
  #Extracting Sequences, identificators, sizes and score.
  VCSeq = list(as.character(datos$aln))
  VCIDs = select(datos$Tabla_unsort, ID, score)
  # This is the size of the read, not the size of the alignment
  SizeList = list(nchar(datos$aln@pattern))
  
  #Bind identificator with sequence and size as to not lose them when we 
  #filter by size
  Graph_Data <- bind_cols(VCIDs, SizeList, VCSeq)
  Graph_Data <- left_join(Graph_Data, datos$Tabla_Original %>% select(ID, Abundance)) %>%
   filter(score >= input$score_threshold) %>% 
   filter(Abundance >= input$abundance_minimum) %>% 
   select(-score) %>% arrange(desc(Abundance))
  
  #Name the tables
  colnames(Graph_Data) = c('ID', 'Size', 'Sequence', 'Abundance')
  
  # Does this make sense?
  #Leave enough margin in order to correct for coverage lax
  Filtered_Graph_Data = filter(Graph_Data, Size < width(datos$Ref)*(1.5))
  Dump = filter(Graph_Data, Size > width(datos$Ref)*(1.5))
  
  TableSize = max(Filtered_Graph_Data['Size'])
  
  datos$unsort_ID_Abundance <- left_join(datos$Tabla_unsort %>% filter(score >= input$score_threshold) %>% select(ID), 
                                         datos$Tabla_raw %>% select(-Freq),  by = 'ID')
  
  #----------------Processing-------------------------
  
  IDsSize = (count(Filtered_Graph_Data))*TableSize
  
  # Data Frame is constructed producing for every ID a variable named
  # Position, Character, and Abundance. Where position is the position
  # of the character, the character represents the base in that
  # position and Abundance the number un sequences with that same base
  # in that position.
  
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
                             select(-Abundance) %>% group_by(Position, Total)
  
  # Commented out since it looks like legacy code to be removed eventually
  
  # LetterPerPosition_Norm = LetterPerPosition_Total %>% 
  #                           ungroup() %>% group_by(Position) %>% 
  #                           mutate_at(vars(Total),funs(./sum(Total))) %>% 
  #                           ungroup(Position) %>% group_by(Chr, Position) %>%
  #                           mutate_at(vars(Total), funs(sum(Total))) %>% 
  #                           group_by(Total, Position) %>% distinct(Chr)
  
  #New and improved LetterPerPosition
  LetterPerPosition_Norm2 = LetterPerPosition_Total %>% 
                             ungroup() %>% group_by(Position) %>% 
                             mutate(Probability = Total/sum(Total)) %>% 
                             ungroup(Position) %>% group_by(Chr, Position) %>%
                             mutate_at(vars(Probability), list(~sum(Probability))) %>% 
                             group_by(Probability, Position) %>% mutate(Total = sum(Total)) %>%  distinct()

  LetterPerPosition_Ref <- tibble(Position = c(1:length(datos$Ref[[1]])))
  LetterPerPosition_Ref$Base <- as.character(as.vector(datos$Ref[[1]]))
  colnames(LetterPerPosition_Ref) <- c('Position','Base')
  
  #For LetterPerPosition_SubChange i want to obtain the accumulated mutation probability for each position.
  LetterPerPosition_SubChange <- left_join(LetterPerPosition_Norm2, LetterPerPosition_Ref)
  
  #With alternative consensus the most abundant mutation for each position is extracted.
  Alternative_Consensus <- LetterPerPosition_SubChange %>% filter(as.character(Chr) != Base) %>% 
                           #filter(as.character(Chr) != '-') %>% 
                           group_by(Position) %>% filter(Probability == max(Probability))
  
  #Deletions_per_position is useful to plot a line graph in conjuntion with BasePerPosition_info 
  Deletion_per_position <- left_join(LetterPerPosition_Ref, LetterPerPosition_SubChange %>% filter(as.character(Chr) == '-') %>% 
                                       select(-Base),by = 'Position') %>% mutate(Chr = '-') %>% 
                                       mutate(Total = case_when(is.na(Total) ~0, TRUE ~as.double(Total))) %>% 
                                       mutate(Probability = case_when(is.na(Probability) ~ 0, TRUE ~ Probability))
  
  #-Insert_per_loci to include insertions in summary plot
  #-Insert locations processing
  
  Insert_data = as.data.frame(indel(datos$aln)@insertion)
  DFGraph_Data <- as.data.frame(Graph_Data)
  Insert_data <- Insert_data %>% 
    mutate(Count = DFGraph_Data['Abundance'][group,]) %>% 
    select(start, end, Count)
  Insert_per_loci <- as.data.frame(c(1:max(width(datos$Sequences))))
  colnames(Insert_per_loci) <- c('Position')
  Insert_per_loci <- Unravel_Positions(Insert_per_loci, Insert_data)
  colnames(Insert_per_loci) <- c('Position', 'Total_Count')
  
  Insert_per_loci <- Insert_per_loci %>% mutate(Probability = Total_Count/sum(datos$Tabla$Abundance))
  #New name to distinguish between purposes (summary plot vs bar chart)
  Insert_per_position <- left_join(LetterPerPosition_Ref, Insert_per_loci, by = 'Position') %>%
                          select(-Total_Count)
  
  
  LetterPerPosition_SubChange <- LetterPerPosition_SubChange %>% filter(as.character(Chr) != '-') %>% 
                                 filter(as.character(Chr) != Base) %>% group_by(Position) %>% mutate(Probability = sum(Probability)) %>% 
                                 select(Position, Base, Probability) %>% distinct()
  
  #Fuse info of LetterPerPosition_SubChange and Alternative_Consensus
  BasePerPosition_info <- left_join(LetterPerPosition_SubChange, Alternative_Consensus %>% select(-Base, -Total, -Probability), by = 'Position')
  # I'm unsure about the purpose of this line, is it to combine positions with different entries?
  BasePerPosition_info <- BasePerPosition_info %>% 
                          mutate(Chr = case_when(Position == lead(Position) ~ paste(Chr, ',' ,lead(Chr), sep = ''), TRUE ~ as.character(Chr))) %>% 
                          distinct(Position, .keep_all = TRUE)
 
  incProgress( 1/8 , detail = 'Deletions')
  
  # Total_Deletions_locations <- Deletion_per_position %>% select(Position, Total) %>% filter(Total != 0)
  
  # Alternative Deletions_locations using alignment data:
  # Deletions_locations <- as_tibble(indel(datos$aln)@deletion) %>% arrange(start) %>% group_by(start) %>% count()
  
  #Plot scale range:
  range = c(0,max(c(Insert_per_loci$Total_Count, Deletion_per_position$Total))+max(c(Insert_per_loci$Total_Count, Deletion_per_position$Total))*0.02)
  
  #-----------------Plotting--------------------------
  
  # plot.new()-error
  
  #-Base Frequency
  
  # colorscheme = c('rgba(255, 77, 77, 0.7)', 'rgba(38, 38, 38, 0.7)', 'rgba(77, 77, 255, 0.7)', 'rgba(128, 255, 128, 0.7)', 'rgba(179, 102, 255, 0.7)')
  # colorscheme = c('red', 'black', 'blue', 'green', 'gray')
  # colorscheme = setNames(colorscheme, c('T', 'G', 'C', 'A', '-'))
  # 
  browser()
  MF <- plot_ly(BasePerPosition_info %>% ungroup(), x = ~Position, y = ~Probability, name = 'Variants',
                type = 'scatter', mode = 'lines', text = ~paste(Chr,'>',Base, ' at Position: ',Position,sep = ''), 
                hovertemplate = paste('<i>Most abundant change</i>: %{text}'), line = list(shape = 'spline', color = "royalblue")) %>% 
        add_trace(xaxis='x2', showlegend = FALSE, opacity = 0) %>%
        add_trace(xaxis='x3', showlegend = FALSE, opacity = 0) %>%
        add_trace(xaxis='x4', showlegend = FALSE, line = list(color = "royalblue"),opacity = 0) %>%
        layout(yaxis = list(range=c(0,1), text='Fraction of variation from reference per position'), 
               xaxis = list(ticktext = BasePerPosition_info$Base[BasePerPosition_info$Base == 'T'], 
                            tickvals = BasePerPosition_info$Position[BasePerPosition_info$Base == 'T'] , 
                            tickmode = "array", 
                            tickangle = 0,
                            title = 'Reference',
                            tickfont=list( color='red' )),
               xaxis2 = list(ticktext = BasePerPosition_info$Base[BasePerPosition_info$Base == 'G'],
                             tickvals = BasePerPosition_info$Position[BasePerPosition_info$Base == 'G'],
                             overlaying='x',
                             tickmode = "array", 
                             tickangle = 0,
                             title = 'Reference',
                             tickfont=list( color='blue' )),
               xaxis3 = list(ticktext = BasePerPosition_info$Base[BasePerPosition_info$Base == 'C'],
                             tickvals = BasePerPosition_info$Position[BasePerPosition_info$Base == 'C'],
                             overlaying='x',
                             tickmode = "array", 
                             tickangle = 0,
                             title = 'Reference',
                             tickfont=list( color='black' )),
               xaxis4 = list(ticktext = BasePerPosition_info$Base[BasePerPosition_info$Base == 'A'],
                             tickvals = BasePerPosition_info$Position[BasePerPosition_info$Base == 'A'],
                             overlaying='x',
                             tickmode = "array", 
                             tickangle = 0,
                             title = 'Reference',
                             tickfont=list( color='green' )))
  
  MF <- MF %>% add_lines(data = Deletion_per_position, x = ~Position, y = ~Probability, name = 'Deletions' ,line = list(shape = 'spline', color = 'tomato'),
                         text = ~paste('del','>',Base, ' at Position: ',Position,sep = ''), hovertemplate = paste('<i>Deletion abundance: </i>%{text}'))
  MF <- MF %>% add_lines(data = Insert_per_position, x = ~Position, y = ~Probability, name = 'Insertions',line = list(shape = 'spline', color = 'limegreen'),
                         text = ~paste('ins','>',Base, ' at Position: ',Position,sep = ''), hovertemplate = paste('<i>Insertion abundance: </i>%{text}'))
  
  output$Mut_Freq <- renderPlotly({print(MF)})
  
  # charts$tmpFileMF <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.pdf'), start = 16), sep = '')
  charts$tmpFileMFhtml <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.html'), start = 16), sep = '')
  charts$tmpFileMFpng <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.png'), start = 16), sep = '')
  tmpFileMFpng <- charts$tmpFileMFpng
  htmlwidgets::saveWidget(MF, paste(CompletePATH,'/', charts$tmpFileMFhtml, sep = ''))
  webshot2::webshot(charts$tmpFileMFhtml, charts$tmpFileMFpng)
  
  #-Deletions per loci
  
  DL <- plot_ly(Deletion_per_position, x = ~Position, y = ~Total, marker = list(color = 'tomato'),
                type = 'bar', name = 'Deletions') %>%
   layout(yaxis = list(title = '# of individual reads', range = range))
  
  output$Del_loc <- renderPlotly({print(DL)})
  
  # charts$tmpFileDL <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.pdf'), start = 16), sep = '')
  # export(DL, charts$tmpFileDL)
  charts$tmpFileDLhtml <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.html'), start = 16), sep = '')
  charts$tmpFileDLpng <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.png'), start = 16), sep = '')
  tmpFileDLpng <- charts$tmpFileDLpng
  htmlwidgets::saveWidget(DL, paste(CompletePATH,'/', charts$tmpFileDLhtml, sep = ''))
  webshot2::webshot(charts$tmpFileDLhtml, charts$tmpFileDLpng)
  
  #-Deletion Sizes
  
  datos$Del_Data = as.data.frame(InDel_Extraction(datos$aln[datos$aln@score >= input$score_threshold], datos$unsort_ID_Abundance)[1])
  
  #Fill data will with ceros for the x-axis to not look too short
  
  rows_to_add = as_tibble()
  if(nrow(datos$Del_Data) < 10){
   rows_to_add <- cbind(as_tibble(c((nrow(datos$Del_Data)+1):10)),as_tibble(rep(0,10-nrow(datos$Del_Data))))
   names(rows_to_add) <- names(datos$Del_Data)
   datos$Del_Data <- bind_rows(datos$Del_Data, rows_to_add)
   tick_distance_Del <- 1
  }else{tick_distance_Del <- 5}
  
  
  DS = plot_ly(datos$Del_Data, y = ~TotalD, x = ~Deletions, type = 'bar', marker = list(color = 'tomato')) %>%
   layout(yaxis = list(title = '# of individual reads'),
          xaxis = list(title = 'Size of deletion extension',
                       dtick = tick_distance_Del,
                       tick0 = 0))
  output$Del_sizes <- renderPlotly({print(DS)})
  
  # charts$tmpFileDS <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.pdf'), start = 16), sep = '')
  # export(DS, charts$tmpFileDS)
  charts$tmpFileDShtml <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.html'), start = 16), sep = '')
  charts$tmpFileDSpng <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.png'), start = 16), sep = '')
  tmpFileDSpng <- charts$tmpFileDSpng
  htmlwidgets::saveWidget(DS, paste(CompletePATH,'/', charts$tmpFileDShtml, sep = ''))
  webshot2::webshot(charts$tmpFileDShtml, charts$tmpFileDSpng)
  
  incProgress( 1/8, detail = 'Insertions')
  
  #-Insertion Sizes
  
  datos$In_Data = as.data.frame(InDel_Extraction(datos$aln[datos$aln@score >= input$score_threshold], datos$unsort_ID_Abundance)[2])
  
  rows_to_add = as_tibble()
  if(nrow(datos$In_Data) < 10){
   rows_to_add <- cbind(as_tibble(c((nrow(datos$In_Data)+1):10)),as_tibble(rep(0,10-nrow(datos$In_Data))))
   names(rows_to_add) <- names(datos$In_Data)
   datos$In_Data <- bind_rows(datos$In_Data, rows_to_add)
   
   if(max(datos$In_Data %>% select(Insertions)) > 10){
    tick_distance_In <- 10
   }else{tick_distance_In <- 1}
  }else{tick_distance_In <- 10}
  
  
  IS = plot_ly(datos$In_Data, y = ~TotalI, x = ~Insertions, type = 'bar', marker = list(color = 'limegreen')) %>%
   layout(yaxis = list(title = '# of individual reads', range = range),
          xaxis = list(title = 'Size of insertion extension',
                       dtick = tick_distance_In,
                       tick0 = 0))
  output$In_sizes <- renderPlotly({print(IS)})
  
  # charts$tmpFileIS <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.pdf'), start = 16), sep = '')
  # export(IS, charts$tmpFileIS)
  charts$tmpFileIShtml <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.html'), start = 16), sep = '')
  charts$tmpFileISpng <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.png'), start = 16), sep = '')
  tmpFileISpng <- charts$tmpFileISpng
  htmlwidgets::saveWidget(IS, paste(CompletePATH,'/', charts$tmpFileIShtml, sep = ''))
  webshot2::webshot(charts$tmpFileIShtml, charts$tmpFileISpng)
  
  #-Insertion per Loci
  incProgress(1/8)
  
  IL = plot_ly(Insert_per_loci, y = ~Total_Count, x = ~Position, type = 'bar', marker = list(color = 'limegreen')) %>%
   layout(yaxis = list(title = '# of individual reads', range = range),
          xaxis = list(range = c(1, 250)))
  output$In_loc <- renderPlotly({print(IL)})
  
  # charts$tmpFileIL <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.pdf'), start = 16), sep = '')
  # export(IL, charts$tmpFileIL)
  charts$tmpFileILhtml <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.html'), start = 16), sep = '')
  charts$tmpFileILpng <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.png'), start = 16), sep = '')
  tmpFileILpng <- charts$tmpFileILpng
  htmlwidgets::saveWidget(IL, paste(CompletePATH,'/', charts$tmpFileILhtml, sep = ''))
  webshot2::webshot(charts$tmpFileILhtml, charts$tmpFileILpng)
  
  output$Dump = renderText(paste('Filter Clusters due to excesive sequence size (>1.5*Reference Sequence length): ',
                                 as.character(count(Dump)),sep = ''))
  # })
 })
 
 if(is.null(datos$Tabla_prev) && nrow(Tabla()) == nrow(datos$Tabla_Original)){
  save(list = c('MF', 'IL','IS','DS','DL', 'Pie_data'), file = paste(session$token,'/default_charts.Rdata', sep = ''))
  save(list = c('tmpFilePiepng', 'tmpFileDSpng', 'tmpFileDLpng',
                'tmpFileISpng', 'tmpFileILpng', 'tmpFileMFpng'), file = paste(session$token,'/defaults_charts_doc.Rdata', sep = ''))
 }
 
 datos$Tabla_prev <- Tabla()
 
})

Edit_Pie_chart <- reactive({
 
 colors <- c('rgb(211,94,96)', 'rgb(128,133,133)', 'rgb(144,103,167)', 
             'rgb(171,104,87)', 'rgb(114,147,203)')
 if(is.null(datos$Pie_data)){return(NULL)}
 
 datos$Pie_data$Groups[-length(datos$Pie_data$Groups)] <- rownames(datos$Tabla)[1:length(datos$Pie_data$Groups)-1]
 
 PC_online <- plot_ly(datos$Pie_data, labels = ~Groups, values = ~Freq, 
                      type = 'pie',
                      textposition = NULL,
                      hoverinfo = 'text',
                      marker = list(colors = colors, 
                                    line = list(color = '#FFFFFF', 
                                                width = 1)), 
                      text = ~Groups, showlegend = FALSE) 
 PC_online <- PC_online %>% layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, 
                                                showticklabels = FALSE),
                                   yaxis = list(showgrid = FALSE, zeroline = FALSE, 
                                                showticklabels = FALSE))})