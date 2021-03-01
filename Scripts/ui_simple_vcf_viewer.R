# Simple vcf viewer ui

# This script includes the user-interface definition of the app.

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