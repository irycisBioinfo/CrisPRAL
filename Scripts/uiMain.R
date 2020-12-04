uiMain <- fluidPage(
 
  # Application Theme
 
 theme = shinytheme("simplex"), 
  useShinyjs(), #To allow JavaScript customization using RCode
  tags$head(tags$script(HTML(press_enter))), #Allow use of enter-key to introduce inputs
# Application title
  titlePanel("All In One"),

  navbarPage( 'Select App', id = 'app', selected = "home",
              ##############################################.
              # Homepage----
              ##############################################.
              tabPanel(
                title = "Home", icon = icon("home"), value = "home",
                mainPanel(width = 9, style="margin-left:4%; margin-right:4%",
                          
                          # Subtitle and help button
                          
                          fluidRow(column(6,(h1("Welcome to the AllInOne tools", style="margin-top:0px;"))), # column end bracket
                                   column(3,actionButton("btn_landing",label="Help: Take tour of the tool",icon=icon('question-circle'),class="down")) # column end bracket
                                   ), # fluidRow end bracket
                          fluidRow(column(3, class="landing-page-column",br(),
                                          actionBttn('go_to_mosaic',
                                            label = 'MOSAIC FINDER',
                                            icon = NULL,
                                            style = "pill",
                                            color = "danger",
                                            size = "lg",
                                            block = TRUE,
                                            no_outline = FALSE)),
                                   column(3, class="landing-page-column",br(),
                                          actionBttn('go_to_mosaic_batch',
                                            label = 'MOSAIC FINDER batch',
                                            icon = NULL,
                                            style = "pill",
                                            color = "warning",
                                            size = "lg",
                                            block = TRUE,
                                            no_outline = FALSE)),
                                   column(3, class="landing-page-column",br(),
                                          actionBttn('go_to_single',
                                            label = 'Single Gene Variant Calling',
                                            icon = NULL,
                                            style = "pill",
                                            color = "primary",
                                            size = "lg",
                                            block = TRUE,
                                            no_outline = FALSE)),
                                          ), # fluidRow end bracket
                ) # mainPanel end bracket
              ), # tabPanel end bracket
             tabPanel( 'Mosaic Finder', value = 'Mosaic_Finder',
                      source('./CrisPRAL_ui.R', local = TRUE)$value ),
             tabPanel( 'Mosaic Finder BATCH', value = 'Mosaic_Finder_batch',
                       source('./ui_Mosaic_batch.R', local = TRUE)$value ),
             tabPanel( 'Lazy Panel Filter', id = 'Lazy_Panel',
                      source('./LazyP_ui.R', local = TRUE)$value ),
             tabPanel('Fastq Quality Analysis', id = 'FastqcRShiny',
                      source('./Fastqc_ui.R', local = TRUE)$value )
  )
)