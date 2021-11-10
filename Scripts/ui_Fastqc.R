#' Run fastQC shiny app
#'
#' @description Returns a shiny app
#' interface to parse many fastQC
#' objects using the package ngsReports
#'
#' @details Currently some plots can
#' take a while to render if the
#' \code{FastqcDataList} passed to
#' \code{fastqcInput} has many elements
#'
#' @param fastqcInput can be a \code{FastqcFileList},
#'  \code{fastqcDataList},
#' or simply a \code{character} vector
#'  of paths to fastqc files.
#'
#'
#' @return UI data for fastQC shiny.
#'
#' @import shinydashboard
#' @import ngsReports
#' @importFrom magrittr %>%
#' @importFrom plotly layout
#' @importFrom plotly plotlyOutput
#' @importFrom plotly renderPlotly
#' @importFrom plotly event_data
#' @importFrom shiny br
#' @importFrom shiny h1
#' @importFrom shiny h5
#' @importFrom shiny observe
#' @importFrom shiny radioButtons
#' @importFrom shiny checkboxInput
#' @importFrom shiny column
#' @importFrom shiny reactiveValues
#' @importFrom shiny observeEvent
#' @importFrom shiny icon
#' @importFrom shiny htmlOutput
#' @importFrom shiny selectInput
#' @importFrom shiny plotOutput
#' @importFrom shiny sliderInput
#' @importFrom shiny renderUI
#' @importFrom shiny renderPrint
#' @importFrom shiny runApp
#' @importFrom shiny textOutput
#' @importFrom shiny renderText
#' @importFrom shiny reactive
#' @importFrom shiny withProgress
#' @importFrom shinyFiles shinyFilesButton
#' @importFrom shinyFiles shinyFileChoose
#' @importFrom shinyFiles shinyDirChoose
#' @importFrom shinyFiles shinyDirButton
#' @importFrom shinyFiles parseFilePaths
#' @importFrom shinyFiles parseDirPath
#' @importFrom shinyFiles getVolumes
#'
#' @examples
#' \dontrun{
#' Get the files included with the package ngsReports
#' packageDir <- system.file("extdata", package = "ngsReports")
#' fileList <- list.files(packageDir, pattern = "fastqc.zip", full.names = TRUE)
#'
#' # Load the FASTQC data as a FastqcDataList object
#' fdl <- getFastqcData(fileList)
#'
#' # Run the Shiny app
#' fastqcShiny(fdl)
#' }
#'
#' @export
#' @rdname fastqcShiny
#'
  
  ## set out menu logic for downstream
  menuItemLogic <- function(flags) {
    ## initiate the menue logic for PASS, WARN, FAIL flags
    ## menuLogic will be a list of length 5 containing information on flags
    menuLogic <- list()
    if (all(flags$Status == "PASS")) {
      menuLogic[[1]] <- "PASS"
      menuLogic[[2]] <- "green"
    }
    else {
      menuLogic[[1]] <- "WARN"
      menuLogic[[2]] <- "yellow"
    }
    
    if (all(flags$Status == "FAIL")) {
      menuLogic[[1]] <- "FAIL"
      menuLogic[[2]] <- "red"
    }
    
    # Fail values
    menuLogic[[3]] <-
      c(sum(flags$Status == "FAIL"), length(flags$Status))
    
    menuLogic[[4]] <-
      c(sum(flags$Status == "WARN"), length(flags$Status))
    
    # Pass values
    menuLogic[[5]] <-
      c(sum(flags$Status == "PASS"), length(flags$Status))
    
    
    
    menuLogic
  }
  
  
  
  ###### renders the box for presence or access of an icon
  renderValBox <- function(count, status, ic, c) {
    renderValueBox({
      valueBox(
        value = paste(count[1], count[2], sep = "/"),
        subtitle = status,
        icon = icon(ic,
                    class = "fa-lg"),
        color = c
      )
    })
  }
  
  
  ### ui page
  dashboardPage(
    dashboardHeader(title = "Fastq Quality Analysis"),
    dashboardSidebar(
      #### gender the menu items
      sidebarMenu(
        menuItem(text = "Summary", tabName = "BS"),
        menuItem(text = "Total Sequences", tabName = "TS"),
        menuItemOutput("BQflag"),
        menuItemOutput("SQflag"),
        menuItemOutput("SCflag"),
        menuItemOutput("GCflag"),
        menuItemOutput("NCflag"),
        menuItemOutput("SLDflag"),
        menuItemOutput("SDLflag"),
        menuItemOutput("OSflag"),
        menuItemOutput("ACflag"),
        menuItemOutput("KCflag")
      ),
      width = 250
    ),
    ## start rendering the dashboard gui
    dashboardBody(
      tabItems(
        tabItem(
          tabName = "BS",
          column(
            width = 2,
            box(h5("Select Input Files:"),
                # shinyFilesButton(
                #  id = "FastqcinputFiles",
                #  label = "Select fastq file",
                #  multiple = FALSE,
                #  title = ''
                fileInput("FastqcinputFiles", "Fastq file", accept = c('.fastq'), multiple = TRUE),
                br(),
                collapsible = TRUE,
                width = NULL,
                title = "Input Files"),
            box(
              h5("Choose FastQC Report:"),
              shinyFilesButton(
                id = "files",
                label = "Select fastq file",
                multiple = TRUE,
                title = ""
              ),
              br(),
              textOutput("report"),
              br(),
              checkboxInput("Sumcluster",
                            "Cluster",
                            value = TRUE),
              collapsible = TRUE,
              width = NULL,
              title = "Options"
            ),
            box(
             downloadButtonModule(id = 'fastqcReport', label = 'Download Fastqc'),
             br(),
             width = NULL,
             title = 'Downloads',
             collapsible = TRUE
            )
          ),
          box(
            h1("Summary of fastQC Flags"),
            h5(
              "Heatmap of fastQC flags
              (pass, warning or fail) for
              each fastQC report"
            ),
            plotlyOutput("SummaryFlags"),
            width = 10
            )
          ),
        tabItem(
          tabName = "TS",
          box(
            checkboxInput("showDup",
                          "Show Duplicated?",
                          value = FALSE),
            collapsible = TRUE,
            width = 2,
            title = "Options"
          ),
          box(
            h1("Total Sequences"),
            h5("Total number of unique and
               duplicated reads in each sample"),
            plotlyOutput("ReadTotals"),
            width = 10
            )
        ),
        tabItem(
          tabName = "BQ",
          column(
            width = 2,
            box(
              radioButtons(
                inputId = "BQplotValue",
                label = "Base Quality",
                choices = c("Mean", "Median"),
                selected = "Mean"
              ),
              checkboxInput("BQcluster",
                            "Cluster",
                            value = TRUE),
              collapsible = TRUE,
              width = NULL,
              title = "Options"
            ),
            valueBoxOutput("BQboxP",
                           width = NULL),
            valueBoxOutput("BQboxW",
                           width = NULL),
            valueBoxOutput("BQboxF",
                           width = NULL)
          ),
          box(
            h1("Per Base Sequence Quality"),
            h5(
              "Per base sequence quality in
              each sample, can either view mean or
              median for each cycle"
            ),
            h5("Click sidebar on heatmap to change
               line plots"),
            plotlyOutput("baseQualHeatmap"),
            br(),
            plotlyOutput("BaseQualitiesSingle"),
            width = 10
            )
          ),
        tabItem(
          tabName = "SQ",
          column(
            width = 2,
            box(
              radioButtons(
                inputId = "SQType",
                label = "Sequence Quality",
                choices = c("Frequency", "Counts"),
                selected = "Frequency"
              ),
              checkboxInput("SQcluster",
                            "Cluster",
                            value = TRUE),
              collapsible = TRUE,
              width = NULL,
              title = "Options"
            ),
            valueBoxOutput("SQboxP",
                           width = NULL),
            valueBoxOutput("SQboxW",
                           width = NULL),
            valueBoxOutput("SQboxF",
                           width = NULL)
          ),
          box(
            h1("Per Sequence
               Quality Scores"),
            h5(
              "Per base sequence quality
              in each sample, can either
              view mean or median for each
              cycle"
            ),
            h5("Click sidebar on heatmap to
               change line plots"),
            plotlyOutput("seqQualHeatmap"),
            br(),
            plotlyOutput("SeqQualitiesSingle"),
            width = 10
            )
            ),
        tabItem(
          tabName = "SC",
          column(
            width = 2,
            box(
              checkboxInput("SCcluster",
                            "Cluster",
                            value = TRUE),
              collapsible = TRUE,
              width = NULL,
              title = "Options"
            ),
            valueBoxOutput("SCboxP",
                           width = NULL),
            valueBoxOutput("SCboxW",
                           width = NULL),
            valueBoxOutput("SCboxF",
                           width = NULL)
          ),
          box(
            h1("Per Base Sequence Content"),
            h5(
              "Per base sequence content in
              each sample, colours at each
              base indicate sequence bias"
            ),
            h5("G = Black, A = Green,
               T = Red, C = Blue"),
            plotlyOutput("SCHeatmap"),
            br(),
            plotlyOutput("SCsingle"),
            width = 10
            )
          ),
        tabItem(
          tabName = "GC",
          column(
            width = 2,
            box(
              checkboxInput("GCcluster",
                            "Cluster",
                            value = TRUE),
              checkboxInput("theoreticalGC",
                            "Normalize To Theoretical GC",
                            value = FALSE),
              htmlOutput("theoreticalGC"),
              htmlOutput("GCspecies"),
              collapsible = TRUE,
              width = NULL,
              title = "Options"
            ),
            valueBoxOutput("GCboxP",
                           width = NULL),
            valueBoxOutput("GCboxW",
                           width = NULL),
            valueBoxOutput("GCboxF",
                           width = NULL)
          ),
          box(
            h1("Per Sequence GC Content"),
            h5(
              "GC content (%) in sample,
              can either view total count
              or frequency"
            ),
            h5("Click sidebar on heatmap
               to change line plots"),
            plotlyOutput("GCheatmap"),
            br(),
            plotlyOutput("GCSingle"),
            width = 10
            )
          ),
        tabItem(
          tabName = "NC",
          column(
            width = 2,
            box(
              checkboxInput("Ncluster",
                            "Cluster",
                            value = TRUE),
              collapsible = TRUE,
              width = NULL,
              title = "Options"
            ),
            valueBoxOutput("NCboxP",
                           width = NULL),
            valueBoxOutput("NCboxW",
                           width = NULL),
            valueBoxOutput("NCboxF",
                           width = NULL)
          ),
          box(
            h1("Per base N content"),
            h5("N content (%) in sample"),
            h5("If dendrogram is truncated
               double click on dendrogram to
               resize"),
            plotlyOutput("NCheatmap"),
            br(),
            plotlyOutput("NCsingle"),
            width = 10
            )
          ),
        tabItem(
          tabName = "SLD",
          column(
            width = 2,
            box(
              radioButtons(
                inputId = "SLType",
                label = "Value to plot",
                choices = c("Frequency",
                            "Counts"),
                selected = "Frequency"
              ),
              checkboxInput("SLcluster",
                            "Cluster",
                            value = TRUE),
              collapsible = TRUE,
              width = NULL,
              title = "Options"
            ),
            valueBoxOutput("SLDboxP",
                           width = NULL),
            valueBoxOutput("SLDboxW",
                           width = NULL),
            valueBoxOutput("SLDboxF",
                           width = NULL)
          ),
          box(
            h1("Sequence Length Distribution"),
            h5(
              "Sequence length distribution
              in each sample, can either
              view total count or frequency"
            ),
            h5("Click sidebar on heatmap to
               change line plots"),
            plotlyOutput("SLHeatmap"),
            br(),
            plotlyOutput("SLSingle"),
            width = 10
            )
          ),
        tabItem(
          tabName = "SDL",
          column(
            width = 2,
            box(
              checkboxInput("Dupcluster",
                            "Cluster",
                            value = TRUE),
              collapsible = TRUE,
              width = NULL,
              title = "Options"
            ),
            valueBoxOutput("SDLboxP",
                           width = NULL),
            valueBoxOutput("SDLboxW",
                           width = NULL),
            valueBoxOutput("SDLboxF",
                           width = NULL)
          ),
          box(
            h1("Sequence Duplication Levels"),
            h5("Sequence duplication in
               each sample"),
            h5("Click sidebar on heatmap
               to change line plots"),
            plotlyOutput("DupHeatmap"),
            br(),
            plotlyOutput("DupSingle"),
            width = 10
            )
          ),
        tabItem(
          tabName = "OS",
          column(
            width = 2,
            box(
              checkboxInput("OScluster",
                            "Cluster",
                            value = TRUE),
              h5("Export Overrepresented
                 Sequences"),
              shinyDirButton(
                id = "dirOS",
                label = "Choose directory",
                title = ""
              ),
              collapsible = TRUE,
              width = NULL,
              title = "Options"
              ),
            valueBoxOutput("OSboxP",
                           width = NULL),
            valueBoxOutput("OSboxW",
                           width = NULL),
            valueBoxOutput("OSboxF",
                           width = NULL)
        ),
        box(
          h1("Overrepresented Sequences"),
          h5("Origin of Overrepresented
             sequences within each sample"),
          plotlyOutput("OSummary"),
          br(),
          plotlyOutput("OSsingle"),
          width = 10
          )
        ),
        tabItem(
          tabName = "AC",
          column(
            width = 2,
            box(
              textInput(
                "ACtype",
                "Regular expression matching the adapter(s) to be plotted",
                value = "Total"
              ),
              checkboxInput("ACcluster",
                            "Cluster",
                            value = TRUE),
              collapsible = TRUE,
              width = NULL,
              title = "Options"
            ),
            valueBoxOutput("ACboxP",
                           width = NULL),
            valueBoxOutput("ACboxW",
                           width = NULL),
            valueBoxOutput("ACboxF",
                           width = NULL)
          ),
          box(
            h1("Adapter content"),
            h5("Adapter content
             (%) across all reads"),
            plotlyOutput("ACheatmap"),
            br(),
            plotlyOutput("ACsingle"),
            width = 10
          )
        ),
        tabItem(
          tabName = "KC",
          column(
            width = 2,
            box(
              checkboxInput("KMcluster",
                            "Cluster",
                            value = TRUE),
              collapsible = TRUE,
              width = NULL,
              title = "Options"
            ),
            valueBoxOutput("KCboxP",
                           width = NULL),
            valueBoxOutput("KCboxW",
                           width = NULL),
            valueBoxOutput("KCboxF",
                           width = NULL)
          ),
          box(
            h1("Kmer Content"),
            h5(
              "Total Identified Kmer
            Count by Position.
            \nPlease load a file to see the top 6 Kmers."
            ),
            plotlyOutput("Kheatmap"),
            br(),
            plotlyOutput("Ksingle"),
            width = 10
          )
        )
  )
  
        )
  )