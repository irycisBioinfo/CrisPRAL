ui <- fluidPage(
  # Dimensions-----
  tags$head(tags$script('
                                var dimension = [0, 0];
                                $(document).on("shiny:connected", function(e) {
                                    dimension[0] = window.innerWidth;
                                    dimension[1] = window.innerHeight;
                                    Shiny.onInputChange("dimension", dimension);
                                });
                                $(window).resize(function(e) {
                                    dimension[0] = window.innerWidth;
                                    dimension[1] = window.innerHeight;
                                    Shiny.onInputChange("dimension", dimension);
                                });
                            ')),
  
  # App title ----
  titlePanel(div("Single Gene Variant Calling", img(src="search.svg")), windowTitle = "Single Gene Variant Calling"), 
  
  # Sidebar panel for inputs ----
  mainPanel(
    wellPanel(
      wellPanel(
        # There's no functionality for single-end right now
        #checkboxInput("single_end", p(strong(" Single end data"))),
        
        h4("Input Files"),
        # fileInput("R1", "Reads 1", accept = c('.fastq', '.fastq.gz')),
        #        conditionalPanel(
        #          condition = "input.single_end == false",
        # fileInput("R2", "Reads 2", accept = c('.fastq', '.fastq.gz'))),
        br(),
        fileInput("dir_VC", "Select .fastq data", accept = c('.gz','fastq'), multiple = TRUE),
        br(),br(),
        fileInput("Reference_VC", "Reference", accept = c('.fasta', '.fastq', '.txt')),
        #        ), # Conditional panel end bracket
        #<-----
      ), # wellPanel closing bracket
      wellPanel(
        h4('Pipeline'),
        selectInput('Aligner', 'Aligner', c('Bowtie2', 'BWA'), selected = 'BWA'),
        selectInput(inputId = 'Variant_Caller', 'Variant Caller', c('GATK', 'freebayes'), selected = 'freebayes'),
        checkboxInput('large_indels', label = ' Look for large indels', value = TRUE),
        br(),
        h4('Options'),
        numericInput('min_len_VC', 'Minimun Read Length (all reads below the set minimun will be removed.)',
                     value = 120, min = 0, max = 900, width = '100%'),
        numericInput('min_qual_VC', label = 'Minimun Quality (minimun average quality of accepted reads in Phred scale.)',
                     value = 30, min = 0, max = 60, step = 10, width = '100%'),
        conditionalPanel( 
          condition = "input.Variant_Caller == 'freebayes'",
        numericInput('min_alternate_fraction', 'Minimun alternate fraction (only for freebayes - minimumn threshold at which allele frequency is considered real)', 
                     value = 0.02, min = 0, max = 0.5,step = 0.1, width = '100%')),
        checkboxInput('remove_duplicates', label = ' Remove duplicates', value = FALSE),
        br(),
        h4('Common Operations'),
        numericInput('pos_correction', 'Translate variant positions (POS column) by a set amount',
                     value = -118, min = -9999, max = 9999, step = 1, width = '100%'),
        checkboxInput('CHROM_transform', label = ' Transform CHROM vcf column to corresponding exons', value = FALSE),
          conditionalPanel( 
            condition = "input.CHROM_transform == true",
            fileInput("exon_mapping", "Exon mapping file", accept = c('.tsv', '.txt')))
      
      ),
      
      wellPanel(actionButton("Accept_VC", "Run"))
    ) # Major wellPanel closing bracket
  ) # Main panel close
) # Fluid Page end bracket