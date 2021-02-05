
#### TODO : MAKE data be saved in specified session #####.

  observeEvent(input$go_to_single, {
    
    updateTabsetPanel(session, 'app', selected = 'Short_Variant_Caller')
    
  })
  
  
  datos <- reactiveValues()
  
  ################################ INPUT DATA ###################################
  ###############################################################################.
  
  Files_location <- reactive({
    #Extract all forward files from the directory specified in input$dir_VC_VC
    
    if(!is.null(input$dir_VC)){
      
      name_datapath <- str_c(dirname(input$dir_VC$datapath),input$dir_VC$name, sep = '/' )
      file.rename(from = input$dir_VC$datapath, to = name_datapath)
      # We need files in the tmp dir to be renamed but nextflow only needs the directory
      return(dirname(input$dir_VC$datapath)[1])
    }else{
      return('Empty')
    }
  })
  
  ################################## PREFIX #####################################
  ###############################################################################.
  
  prefix_VC <- reactive({
    
    # Indentify prefixes for sample identification
    
    if(!is.null(R1_file())){ 
      
      #Split the path from input$dir_VC/paths to obtain only the last part of the path.
      
      return(str_split(R1_file()$name, '\\_', simplify = TRUE)[,1])
    }else{return(NULL)}
    
  })  
  
  ################################## Exon Mapper ###############################
  ##############################################################################.
  ## TODO : provide an example of the mapping file structure  ###
  
  exon_mapping <- reactive({
    
    if(input$CHROM_transform == TRUE){
      return(paste('--exon_annotation',input$exon_mapping$datapath))
    }else{return('')}
    
  })
  
  ################################## Min Alternate Fraction ####################
  ##############################################################################.
  
  min_alternate_fraction <- reactive({
    if(input$Variant_Caller == 'freebayes'){
    return(str_c('--min_alt_fraction', input$min_alternate_fraction))
    }else{return('')}
    
  })
  
  ################################## POS correction ############################
  ##############################################################################.
  
  pos_correction <- reactive({
    
      return(paste('--position_correction',input$pos_correction))
    
  })
  
  ################################## Large indels flag ##########################
  ##############################################################################.
  
  
  large_indels <- reactive({input$large_indels})
  
  ################################## RUN #######################################
  ##############################################################################.
  
  observeEvent(input$Accept_VC,{
    
    req(input$dir_VC, input$Reference_VC)
    
    #Checks for required files and displays a dialog warning its missing.
    #Makes no sense to check if the required values are present inside a function that will only execute if the required values are present.
    required_values <- list("A Reference sequence is" = input$Reference_VC$name, 'Fastq files are' = input$dir_VC)
    name <- names(required_values[required_values == 'NULL'])
    datos$tmppipelinedir <- str_c('data', str_sub(tempfile(), start = 21))
    dir.create(datos$tmppipelinedir)
    
    if( any(required_values == 'NULL') ) {
      
      name <- names(required_values[required_values == 'NULL'])
      
      showModal(modalDialog(
        fluidPage(
          h1(strong("Warning"),align="center"),
          hr(),
          fluidRow(column(4,offset = 4,
                          div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                              img(src="caution-icon.png", height='150', width='150', align="middle")))),
          h2(paste(name ," required to align reads." ,sep = ''), align = 'center'),
          easyClose = TRUE,
          footer = NULL
        )))
      
    }
    
    command_aligner = paste('nextflow',
                            './Scripts/nextflow_GenomeMapper.nf',
                            '--indir',
                            Files_location(),
                            '--genome',
                            input$Reference_VC$datapath,
                            '--adapters',
                            '/home/bioinfo/Proyecto_NF1/adapters/TruSeq3-PE-2.fa',
                            '--paired',
                            '--aln',
                            str_to_lower(input$Aligner),
                            '--remove_duplicates',
                            input$remove_duplicates,
                            '--threads 4',
                            '--outdir',
                            str_c('./',as.character(datos$tmppipelinedir))
    )
    
    command_SNP = paste('nextflow',
                        './Scripts/nextflow_jointVariantCalling_NF1.nf',
                        '--GVCFmode',
                        'false',
                        '--indir',
                        Files_location(),
                        '--genome',
                        input$Reference_VC$datapath,
                        '--vc',
                        input$Variant_Caller,
                        min_alternate_fraction(),
                        exon_mapping(),
                        pos_correction(),
                        '--threads 4',
                        '--outdir',
                        str_c(as.character(datos$tmppipelinedir))
    )
    
    command_INDEL = paste('nextflow',
                          './Scripts/nextflow_LargeIndels.nf',
                          '--indir',
                          Files_location(),
                          '--genome',
                          input$Reference_VC$datapath,
                          '--adapters',
                          '/home/bioinfo/Proyecto_NF1/adapters/TruSeq3-PE-2.fa',
                          '--paired',
                          '--remove_duplicates',
                          input$remove_duplicates,
                          min_alternate_fraction(),
                          pos_correction(),
                          exon_mapping(),
                          '--threads 4',
                          '--outdir',
                          str_c(as.character(datos$tmppipelinedir))
    )

    system(command_aligner)
    system(command_SNP)
    if(large_indels()){
      system(command_INDEL)
    }
    #}
    
    setwd(datos$tmppipelinedir)
    final_files <- c(dir('.', pattern = 'final_variant_calling_files', full.names = TRUE),
                     dir('.', pattern = "Large_indels_variant_calling_files", full.names = TRUE))
    
    if(input$remove_duplicates){
      final_files <- c(final_files, dir('alignment', pattern = '*dups.sort.bam', full.names = TRUE))
    }else{final_files <- c(final_files, dir('alignment', pattern = '*bowtie2.sort.bam', full.names = TRUE))}
    
    path_to_zip_file <- str_c(CompletePATH,'/',as.character(datos$tmppipelinedir),'/','results.zip')
    incProgress( 1/4 ,detail = 'zipping files...')
    zip(path_to_zip_file, files = final_files)
    setwd(CompletePATH)
    
    showModal(modalDialog(
      fluidPage(
        h1(strong("File ready"),align="center"),
        hr(),
        fluidRow(column(4,offset = 4,
                        div(style='max-width:150px;max-height:150px;width:3%;height:5%;', 
                            img(src='ready_for_download.png', height='150', width='150', align="middle")))),
        h2(paste("Click to download file." ,sep = ''), align = 'center'),
        fluidRow(column(4,offset = 4.5,
                        callModule(downloadZIP, 'downloadZIP_NF1', data = path_to_zip_file, 'results'))),
        easyClose = TRUE,
        footer = NULL
      )))
    # 
    # output$vcf <- renderDT({datos$vcf_df})
    # output$vcftable <- renderUI({DTOutput('vcf', width = str_c(input$dimension[1]/4, 'px', sep = ''))})
    # 
  }) # observeEvent RUN end brackets