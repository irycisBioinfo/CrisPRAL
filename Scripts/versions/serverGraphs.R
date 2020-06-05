#Graph generating script.

observeEvent(input$GG,{
 
 withProgress(message = "Plotting graphs: ",
              detail = 'Calculating', {
               #-Pie-Chart-processing----
               
               Other <- datos$Tabla %>% filter(Freq < 1)
               datos$Pie_data <- select(datos$Tabla, ID, Freq) %>% filter(Freq > 1) %>% 
                add_row(ID = paste('<1% groups', 
                                   '(', count(Other['Freq']) , ')', 
                                   sep = ''), 
                        Freq = sum(Other['Freq']))
               
               #-Pie-Chart-plotting online
               
               Groups <- Group_list(datos$Pie_data) #- Function in Functions.R module
               datos$Pie_data <- add_column(datos$Pie_data, Groups)
               
               output$Pie_summary <- renderText({print(Edit_Pie_chart())})
               
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
               
               charts$tmpFilePiehtml <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.html'), start = 16), sep = '')
               charts$tmpFilePiepng <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.png'), start = 16), sep = '')
               htmlwidgets::saveWidget(PC_pdf, paste(CompletePATH,'/', charts$tmpFilePiehtml, sep = ''))
               webshot::webshot(charts$tmpFilePiehtml, charts$tmpFilePiepng)
               
               output$Pie_summary <- renderPlotly({Edit_Pie_chart()})
               
               #-----------Setting up the Data for further plotting----------------
               
               incProgress( 1/8 ,detail = 'Mutation Frequency') #User feedback
               
               #Extracting Sequences, identificators, sizes
               VCAbundance = select(datos$Tabla_Original, Abundance, score)
               VCSeq = list(as.character(datos$aln@pattern))
               VCIDs = select(datos$Tabla_Original, ID)
               SizeList = list(nchar(datos$aln@pattern))
               
               #Bind identificator with sequence and size as to not lose them when we 
               #filter by size
               Graph_Data = bind_cols(VCIDs, SizeList, VCSeq, VCAbundance) %>% 
                filter(score >= input$score_threshold) %>% select(-score)
               
               #Name the tables
               colnames(Graph_Data) = c('ID', 'Size', 'Sequence', 'Abundance')
               
               #Leave enough margin in order to correct for coverage lax
               Filtered_Graph_Data = filter(Graph_Data, Size < 250*(2-input$cov))
               Dump = filter(Graph_Data, Size > 250*(2-input$cov))
               
               TableSize = max(Filtered_Graph_Data['Size'])
               
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
               
               LetterPerPosition_Norm = LetterPerPosition_Total %>% 
                ungroup() %>% group_by(Position) %>% 
                mutate_at(vars(Total),funs(./sum(Total))) %>% 
                ungroup(Position) %>% group_by(Chr, Position) %>%
                mutate_at(vars(Total), funs(sum(Total))) %>% 
                group_by(Total, Position) %>% distinct(Chr)
               
               #-Isolating deletions
               
               incProgress( 1/8 ,detail = 'Deletions')
               
               Deletions_locations <- LetterPerPosition_Total %>% filter(Chr == '-') %>%
                group_by(Position) %>%
                mutate(Total = sum(Total)) %>% ungroup() %>%
                select(Position, Total) %>% distinct()
               
               
               # Alternative Deletions_locations using alignment data:
               # Deletions_locations <- as_tibble(indel(datos$aln)@deletion) %>% arrange(start) %>% group_by(start) %>% count()
               
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
               
               #-----------------Plotting--------------------------
               
               # plot.new()-error
               
               #-Base Frequency
               
               colorscheme = c('red', 'black', 'blue', 'green', 'purple')
               colorscheme = setNames(colorscheme, c('T', 'G', 'C', 'A', '-'))
               
               #Determinaing the range of xaxis en MF plot
               min_value = 1
               max_value = width(datos$Ref@ranges)
               
               MF <- plot_ly(LetterPerPosition_Norm, x = ~Position, y = ~Total, 
                             type = 'bar', name = ~Chr, color= ~Chr , 
                             colors = colorscheme) %>%
                layout(xaxis = list(range = range(min_value, max_value)), yaxis = list(title = 'Count'), barmode = 'stack')
               output$Mut_Freq <- renderPlotly({print(MF)})
               
               # charts$tmpFileMF <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.pdf'), start = 16), sep = '')
               charts$tmpFileMFhtml <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.html'), start = 16), sep = '')
               charts$tmpFileMFpng <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.png'), start = 16), sep = '')
               # export(MF, charts$tmpFileMF)
               
               htmlwidgets::saveWidget(MF, paste(CompletePATH,'/', charts$tmpFileMFhtml, sep = ''))
               webshot::webshot(charts$tmpFileMFhtml, charts$tmpFileMFpng)
               
               #-Deletions per loci
               
               DL <- plot_ly(Deletions_locations, x = ~Position, y = ~Total, 
                             type = 'bar',
                             name = 'Deletions') %>%
                layout(yaxis = list(title = 'Count'))
               output$Del_loc <- renderPlotly({print(DL)})
               
               # charts$tmpFileDL <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.pdf'), start = 16), sep = '')
               # export(DL, charts$tmpFileDL)
               charts$tmpFileDLhtml <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.html'), start = 16), sep = '')
               charts$tmpFileDLpng <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.png'), start = 16), sep = '')
               
               htmlwidgets::saveWidget(DL, paste(CompletePATH,'/', charts$tmpFileDLhtml, sep = ''))
               webshot::webshot(charts$tmpFileDLhtml, charts$tmpFileDLpng)
               
               #-Deletion Sizes
               
               datos$Del_Data = as.data.frame(InDel_Extraction(datos$aln, datos$unsort_ID_Abundance)[1])
               
               DS = plot_ly(datos$Del_Data, y = ~TotalD, x = ~Deletions, type = 'bar' ) %>%
                layout(yaxis = list(title = 'Count'),
                       xaxis = list(title = 'Size'))
               output$Del_sizes <- renderPlotly({print(DS)})
               
               # charts$tmpFileDS <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.pdf'), start = 16), sep = '')
               # export(DS, charts$tmpFileDS)
               charts$tmpFileDShtml <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.html'), start = 16), sep = '')
               charts$tmpFileDSpng <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.png'), start = 16), sep = '')
               
               htmlwidgets::saveWidget(DS, paste(CompletePATH,'/', charts$tmpFileDShtml, sep = ''))
               webshot::webshot(charts$tmpFileDShtml, charts$tmpFileDSpng)
               
               incProgress( 1/8 ,detail = 'Inserts')
               
               #-Insertion Sizes
               
               datos$In_Data = as.data.frame(InDel_Extraction(datos$aln, datos$unsort_ID_Abundance)[2])
               
               IS = plot_ly(datos$In_Data, y = ~TotalI, x = ~Insertions, type = 'bar', color = 'rgba(255,0,0,1)') %>%
                layout(yaxis = list(title = 'Amount of Reads'),
                       xaxis = list(title = 'Size'))
               output$In_sizes <- renderPlotly({print(IS)})
               
               # charts$tmpFileIS <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.pdf'), start = 16), sep = '')
               # export(IS, charts$tmpFileIS)
               charts$tmpFileIShtml <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.html'), start = 16), sep = '')
               charts$tmpFileISpng <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.png'), start = 16), sep = '')
               
               htmlwidgets::saveWidget(IS, paste(CompletePATH,'/', charts$tmpFileIShtml, sep = ''))
               webshot::webshot(charts$tmpFileIShtml, charts$tmpFileISpng)
               
               #-Insertion per Loci
               incProgress(1/8)
               
               IL = plot_ly(Insert_per_loci, y = ~Total_Count, x = ~Position, 
                            type = 'bar', color = 'rgba(255,0,0,1)') %>%
                layout(yaxis = list(title = 'Amount of Reads'),
                       xaxis = list(range = c(1, 250)))
               output$In_loc <- renderPlotly({print(IL)})
               
               # charts$tmpFileIL <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.pdf'), start = 16), sep = '')
               # export(IL, charts$tmpFileIL)
               charts$tmpFileILhtml <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.html'), start = 16), sep = '')
               charts$tmpFileILpng <- paste(datos$tmppipelinedir, str_sub(tempfile(fileext = '.png'), start = 16), sep = '')
               
               htmlwidgets::saveWidget(IL, paste(CompletePATH,'/', charts$tmpFileILhtml, sep = ''))
               webshot::webshot(charts$tmpFileILhtml, charts$tmpFileILpng)
               
               #- Dumped count
               
               output$Dump = renderText(paste('Filter Clusters due to sequence size limits: ',
                                              as.character(count(Dump)),sep = ''))
               # })
              })})

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