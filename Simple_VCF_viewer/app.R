#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyBS)
library(vcfR)
library(tidyverse)
library(DT)

ScriptsPATH = paste(getwd(),'/..', '/Scripts', sep = '')
source( str_c( ScriptsPATH, '/Modules/downloadFileModule.R'))
source( str_c( ScriptsPATH, '/Modules/Functions.R' ))

options( shiny.maxRequestSize = 500 * 1024 ^ 2 )

ui <- fluidPage(

    # Application title
    titlePanel("Simple VCF viewer"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            wellPanel(
                # checkboxInput(inputId = 'multi_sample', label = tags$s('Check for vcf with multiple samples'), value = FALSE),
                
                uiOutput(outputId = 'input_vcf_variable'),
                bsPopover('input_vcf_variable', title = "Input Data" ,
                          content = 'Select all your fastq or fastq.gz files, make sure pairs are defined with the basic nomenclature in the file names: _R1_ & _R2_', 
                          placement = "bottom", trigger = "hover", options = NULL)
                
            ),  # Well Panel end bracket
            wellPanel(
                
                popify(fileInput('vcf_reference', label = 'Input reference sequence', 
                                 multiple = FALSE, accept = c('fa','fasta','txt'), placeholder = 'Select file'),
                       "Input reference sequence",'For the moment reference is only used to determine the variant names, is does not need to be the one used for the variant calling, though it probably makes sense that it is for the variant position to make sense .'),
                br(),
                actionButton(inputId = 'run_vcf_reader', label = 'Load')
                
            )
        ),  # SidebarPanel end bracket

        # Show a plot of the generated distribution
        mainPanel(
            h3('Variant Calling File:'),
            DTOutput(outputId = 'vcf_table'),
            verbatimTextOutput('hoverIndex'),
            h3('Common Operations:',icon('tools')),
            fluidRow(
                column(width = 3,
                       wellPanel(
                           
                           h4('Download selected variants'),
                           checkboxInput(inputId = 'select_all_variants', label = 'Select all', value = FALSE),
                           downloadButtonModule("download_variants", "Download"),
                           
                       ) # Well Panel end bracket
                       ),# Column end bracket
                column(width = 6,
                       wellPanel(
                           
                            h4('Selected Variant:'),
                            textOutput(outputId = 'display_variant', inline = TRUE)
                           
                       )) # Well Panel end bracket
                       ),
            fluidRow(
                column(width = 5,
                    wellPanel(
                        h4('Recalculate Alternate Frequency (AF) values',icon('question-circle'),id='Recalc_AF'),
                        bsPopover('Recalc_AF', title = "Recalculate Alternate Frequency (AF) values" ,content = 'Alternate frequency commonly needs to be recalculated, some variant callers assume diploid or haploid distribution of the data, and in consequence they print 0,0.5 or 1 in the AF column. Despite the real AF value being something in between.',
                                  placement = "top", trigger = "hover", options = NULL),
                        br(),
                        actionButton('Recalculate_AF',label = 'Recalculate AF', icon = icon("redo-alt"))
                )),# wellPanel & column

                column(width = 5,
                   wellPanel(
                       h4('Fixed translation',icon('question-circle'),id='Fixed_trans'),
                       bsPopover('Fixed_trans', title = "Fixed translation" ,content = 'Translate the selected columns by a fixed amount, useful to map vcf to a different reference.',
                                 placement = "top", trigger = "hover", options = NULL),
                       numericInput('translate_int', label = '', value = 0,min = -9999999, max = 9999999),
                       br(),
                       actionButton('Translate',label = 'Translate', icon = icon("cut"))
                   )),# wellPanel & column

                column(width = 5,
                       wellPanel(
                           h4('Map position',icon('question-circle'),id='pos_map_info'),
                           bsPopover('pos_map_info', title = "Map position" ,content = HTML('You can annotate the vcf using a two column table that maps a range of position to a characteristic, such as an exon. The file MUST contain an ANNOTATION, a START and an END column. Position ranges (START:END) must not overlap and each range must map to a single ANNOTATION entry.'),
                                     placement = "top", trigger = "hover", options = NULL),
                           fileInput("pos_mapping", "", accept = c('.tsv', '.txt')),
                           # br(),
                           actionButton('Go_mapping',label = 'Go', icon = icon("random"))
                       ))# wellPanel & column
                  )# fluidRow  end bracket
            

        ) # Main Panel end bracket
    ) # SideBar layout end bracket
) # fluidPage end bracket

# Define server logic required to draw a histogram
server <- function(input, output) {
    datos <- reactiveValues()
    
    output$infocircle <- renderImage({outfile <- "ready_for_download.png"
        list(src = outfile,
            contentType = 'image/png',
            width = 400,
            height = 300,
            alt = "This is alternate text")
    }, deleteFile=FALSE)
    
    observeEvent(input$run_vcf_reader, {
        
        vcf_r <- reactive({read.vcfR(input$vcf_file$datapath)})
        
        datos$vcf_tidy <- vcfR2tidy(vcf_r())
        # Fix numeric interpretation of data:
        datos$vcf_tidy$fix <- datos$vcf_tidy$fix %>% mutate(across(.cols = c(AC,AO,AF,PAO,QA,PQA,SAF,RPP,RPL,RPR,EPP,DPRA,MEANALT,MQM,PAIRED,PAIREDR), as.numeric))
        
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
        
        info_columns <- which(colnames(datos$vcf_tidy$fix) %in% datos$vcf_tidy$meta$ID)
        basic_columns <- which(! colnames(datos$vcf_tidy$fix) %in% datos$vcf_tidy$meta$ID)
        
        output$vcf_table <- DT::renderDT(extensions = c('Buttons','FixedColumns','Scroller'),
                                     filter = 'top',
                                     options = list(dom = 'Bfrtip',
                                                    scrollX = TRUE, fixedColumns = list(leftColumns = 2),
                                                    deferRender = TRUE, scrollY = 400, scroller = TRUE,
                                                    buttons = list('copy', list(extend = 'collection',
                                                                               buttons = c('csv', 'excel'),
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
        
        datos$reference_type_ab <- switch(reference_type(), 'NF' = 'c', 'NC' = 'g', 'NP' = 'p', 'NG' = 'g', 'NR' = 'c', 'LR' = 'g')
        
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

############# VARIANTS NOMENCLATURE ###############
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
    
    observeEvent(input$Go_mapping,{
        
        if(isFALSE(check_annot_data())){return()}
        
        vcf_pos <- datos$vcf_tidy$fix$POS
        vcf_annot <- c()
        for(entry in vcf_pos){
            annot.idx <- dplyr::last(which(entry > gene_annotation()$START))
            annot <- gene_annotation()$ANNOTATION[annot.idx]
            vcf_annot <- c(vcf_annot, annot)
        }
        
        datos$vcf_tidy$fix <- bind_cols(datos$vcf_tidy$fix, POS_ANNOTATION = vcf_annot) %>% relocate(POS_ANNOTATION,.after = POS)
        
        })
    
} # server end bracket


# Run the application 
shinyApp(ui = ui, server = server)
