# Isolated script for the loading of the simple_vcf_viewer operations

  observeEvent(input$go_to_simple_vcf, {
    
    updateTabsetPanel(session, 'app', selected = 'Simple_VCF')
    
  })
  
  output$infocircle <- renderImage({outfile <- "ready_for_download.png"
  list(src = outfile,
       contentType = 'image/png',
       width = 400,
       height = 300,
       alt = "This is alternate text")
  }, deleteFile=FALSE)
  
  ################ RUN ####################
  ########################################.
  
  observeEvent(input$run_vcf_reader, {
    
    vcf_r <- reactive({read.vcfR(input$vcf_file$datapath)})
    
    datos$vcf_tidy <- vcfR2tidy(vcf_r())
    
    ####
    # Fix numeric interpretation of data:
    numeric_INFO_values <- datos$vcf_tidy$meta %>% filter(Tag == 'INFO') %>% filter(Type %in% c('Integer', 'Float')) %>% pull(ID)
    datos$vcf_tidy$fix <- datos$vcf_tidy$fix %>% mutate(across(.cols = numeric_INFO_values , as.numeric))
    ####
    
    vcf.meta <- reactive({datos$vcf_tidy$meta})
    vcf.table <- reactive({datos$vcf_tidy$fix})
    vcf.gt <- reactive({datos$vcf_tidy$gt})
    
    vcf_reference <- reactive({
      if(is.null(input$vcf_reference)){return('no-reference')
      }else{return(input$vcf_reference$name)}
    })
    
    reference_type <- reactive({return(str_sub(vcf_reference(), 1,2))})
    reference_ID <- reactive({return(str_split(vcf_reference(), '\\.')[[1]][1])})
    
    output$hoverIndex <- renderText({
      # Find column user is hovering (we add one because on display we are hiding column #1)
      HoverColumnName <- colnames(datos$vcf_tidy$fix)[input$hoverIndexJS+1]
      
      if(is_empty(HoverColumnName)){UI_out <- paste('Hover over an INFO column to display column description')}else{
        
        # Find position in "meta" table
        INFO_def_loc <- match(HoverColumnName, datos$vcf_tidy$meta$ID)
        # Query meta table for definition:
        if(!is.na(INFO_def_loc)){
          
          UI_out <- paste(datos$vcf_tidy$meta$ID[INFO_def_loc],'-',datos$vcf_tidy$meta$Description[INFO_def_loc], sep = ' ')
          
        }else{UI_out <- paste('Hover over an INFO or column to display column description')}}
      return(paste("INFO column details:", UI_out))
      
      # return(paste("displaying:",input$hoverIndexJS))
    })
    
    info_columns <- which(colnames(datos$vcf_tidy$fix %>% select(-ChromKey)) %in% datos$vcf_tidy$meta$ID)
    basic_columns <- which(! colnames(datos$vcf_tidy$fix %>% select(-ChromKey)) %in% datos$vcf_tidy$meta$ID)
    
    #### LOAD VCF TABLE ####
    
    output$vcf_table <- DT::renderDT(extensions = c('Buttons','FixedColumns','Scroller'),
                                     filter = 'top',
                                     options = list(dom = 'Bfrtip',
                                                    scrollX = TRUE, fixedColumns = list(leftColumns = 2),
                                                    deferRender = TRUE, scrollY = 400, scroller = TRUE,
                                                    buttons = list('copy', list(extend = 'collection',
                                                                                buttons = list(
                                                                                  list(extend = 'csv', filename = 'VCF_file'),
                                                                                  list(extend = 'excel', filename = 'VCF_file')),
                                                                                text = 'Download'),
                                                                   'colvis', list(extend = 'colvisGroup',
                                                                                  text = 'INFO',
                                                                                  show = basic_columns,
                                                                                  hide = info_columns),
                                                                   list(extend = 'colvisGroup',
                                                                        text = 'DEFAULT',
                                                                        show = info_columns))),
                                     selection = list('multiple',target = 'row+column'),
                                     callback = JS("
                                                /* code for columns on hover */
                                                table.on('mouseenter', 'td', function() {
                                                    var td = $(this);
                                                    var info_out = table.cell( this ).index().column;
                                                    Shiny.onInputChange('hoverIndexJS', info_out);
                                                });"),{
                                                  vcf.table() %>% select(-ChromKey)
                                                }) #output end brackets
    
    datos$reference_type_ab <- switch(reference_type(), 'NF' = 'c', 'NC' = 'g', 'NP' = 'p', 'NG' = 'g', 'NR' = 'c', 'LR' = 'g', 'NM' = 'r')
    
    if(!is.null(datos$reference_type_ab)){
      output$display_variant <- renderText(str_c(reference_ID(),':',
                                                 datos$reference_type_ab,'.',
                                                 vcf.table()[input$vcf_table_row_last_clicked,]$POS,
                                                 vcf.table()[input$vcf_table_row_last_clicked,]$REF,'>',
                                                 vcf.table()[input$vcf_table_row_last_clicked,]$ALT))
    }else{
      output$display_variant <- renderText(str_c(vcf.table()[input$vcf_table_row_last_clicked,]$POS,
                                                 vcf.table()[input$vcf_table_row_last_clicked,]$REF,
                                                 '>',vcf.table()[input$vcf_table_row_last_clicked,]$ALT))
    }
  })#observe event end brackets
  
############# TABLE PROXY #######################
  #################################.
  
  
  vcf_table_proxy <- DT::dataTableProxy("vcf_table")
  observeEvent(input$select_all_variants, {
    
    if (isTRUE(input$select_all_variants)) {
      DT::selectRows(vcf_table_proxy, input$vcf_table_rows_all)
    } else {
      DT::selectRows(vcf_table_proxy, NULL)
    }
    
  }) # observeEvent end bracket
  
############ VARIANTS NOMENCLATURE ###############
  #################################.
  
  variants_to_download <- reactive({
    if(!is.null(datos$reference_type_ab)){
      data <- tibble(str_c(reference_ID(),':',
                           datos$reference_type_ab,'.',
                           datos$vcf_tidy$fix[input$vcf_table_rows_selected,]$POS,
                           datos$vcf_tidy$fix[input$vcf_table_rows_selected,]$REF,'>',
                           datos$vcf_tidy$fix[input$vcf_table_rows_selected,]$ALT)) %>% pull()}
    else{
      data <- tibble(str_c(
        datos$vcf_tidy$fix[input$vcf_table_rows_selected,]$POS,
        datos$vcf_tidy$fix[input$vcf_table_rows_selected,]$REF,'>',
        datos$vcf_tidy$fix[input$vcf_table_rows_selected,]$ALT)) %>% pull()
    }
    
    return(data)
    
  })
  
  ############# INPUT VCF UI #######################
  #################################.
  
  output$input_vcf_variable <- renderUI({
    # Change multiple = FALSE for multiple = input$multi_sample when ready.
    fileInput(inputId = 'vcf_file', label = 'Input vcf', multiple = FALSE, accept = c('.vcf','.tsv'), placeholder = 'Select vcf file')
  })
  
  ############# DOWNLOAD FILES #######################
  #################################.
  variant_file_name <- reactive({'variant_list.txt'})
  callModule(downloadFile, 'download_variants', data = variants_to_download, name = variant_file_name)
  
  ############# RECALCULATE AF #######################
  #################################.
  observeEvent(input$Recalculate_AF, {datos$vcf_tidy$fix <- datos$vcf_tidy$fix %>% mutate(AF = signif(AO/DP,digits = 2))})# observeEvent end bracket
  
  ############# FIXED TRANSLATION ####################
  ################################.
  
  observeEvent(input$Translate, {
    
    if(is.null(input$vcf_table_columns_selected)){
      uiOutput( 
        showModal(modalDialog(
          
          fluidPage(
            h1(strong("Warning"),align="center"),
            hr(),
            fluidRow(column(4,offset = 4,
                            div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                                img(src="caution-icon.png", height='100', width='100', align="middle")))),
            h2("You must at least select one column.", align = 'center'),
            easyClose = TRUE,
            footer = NULL
          )))) # uiOutput end bracket
    }else if(is.numeric(datos$vcf_tidy$fix[,input$vcf_table_columns_selected+1] %>% pull())){
      
      # replaceData(vcf_table_proxy, datos$vcf_tidy$fix %>% mutate(across(.cols = c(input$vcf_table_columns_selected)+1, ~.+input$translate_int)), clearSelection = 'all')
      datos$vcf_tidy$fix <- datos$vcf_tidy$fix %>% mutate(across(.cols = c(input$vcf_table_columns_selected)+1, ~.+input$translate_int))
      
    }
    # Displays a warning in case a non-numeric column is selected.
    else{uiOutput( 
      showModal(modalDialog(
        
        fluidPage(
          h1(strong("Warning"),align="center"),
          hr(),
          fluidRow(column(4,offset = 4,
                          div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                              img(src="caution-icon.png", height='100', width='100', align="middle")))),
          h2("Only select numeric columns.", align = 'center'),
          easyClose = TRUE,
          footer = NULL
        ))) # modal dialog end bracket
    ) # uiOutput end bracket
    }
  })# observeEvent end bracket
  
############# POSITION MAPING ######################
  ################################.
  
  #### ASSERTS ####
  
  gene_annotation <- reactive({
    read_tsv(input$pos_mapping$datapath, col_names = TRUE)
  })
  
  check_annot_data <- eventReactive(gene_annotation(),{
    
    mandatory_columns = c('ANNOTATION', 'START', 'END')
    if(sum(mandatory_columns %in% colnames(gene_annotation())) != 3){
      
      showModal(modalDialog(
        fluidPage(
          h1(strong("Warning"),align="center"),
          hr(),
          fluidRow(column(4,offset = 4,
                          div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                              img(src="caution-icon.png", height='100', width='100', align="middle")))),
          h2(paste(paste(mandatory_columns[!mandatory_columns %in% colnames(gene_annotation())], collapse = ' and '),"could be missing"), align = 'center'),
          easyClose = TRUE,
          footer = NULL
        )))# modal dialog end bracket
      return(FALSE)
    }
    
    # Assert position ranges DO NOT overlap
    if(isTRUE(any(lag(gene_annotation()$END) > gene_annotation()$START))){
      showModal(modalDialog(
        fluidPage(
          h1(strong("Warning"),align="center"),
          hr(),
          fluidRow(column(4,offset = 4,
                          div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                              img(src="caution-icon.png", height='100', width='100', align="middle")))),
          h2("Position ranges cannot overlap. i.e map to multiple annotations.", align = 'center'),
          easyClose = TRUE,
          footer = NULL
        )))# modal dialog end bracket
      return(FALSE)
    }
    return(TRUE)
  })
  
  #### EVENT ####
  
  observeEvent(input$Go_mapping,{
    
    if(isFALSE(check_annot_data())){return()}
    # check if pos_annotation already exists
    if(any(str_detect(colnames(datos$vcf_tidy$fix), pattern = 'POS_ANNOTATION'))){
      datos$vcf_tidy$fix <-  datos$vcf_tidy$fix %>% select(-POS_ANNOTATION)
          }
    
    vcf_pos <- datos$vcf_tidy$fix$POS
    vcf_annot <- c()
    for(entry in vcf_pos){
      annot.idx <- dplyr::last(which(entry > gene_annotation()$START))
      annot <- gene_annotation()$ANNOTATION[annot.idx]
      vcf_annot <- c(vcf_annot, annot)
    }
    
    datos$vcf_tidy$fix <- bind_cols(datos$vcf_tidy$fix, POS_ANNOTATION = vcf_annot) %>% relocate(POS_ANNOTATION,.after = POS)
    
  })
