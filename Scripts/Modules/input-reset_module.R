inputReset_UI <- function(id, label){
  
  ns <- NS(id)
  tagList(
    useShinyjs(),
    fileInput(ns('file'), label, accept = c('.fasta','.txt')),
    actionButton(ns('resetButton'), 'Reset')
         )
  }


inputReset_server <- function(input, output, session){
  # ns <- session$ns
  
  file_upload <- reactiveValues(upload_Target_seq = NULL)
  
  observeEvent(input$file, {
    file_upload$file_upload_sate <- 'uploaded'})
  observeEvent(input$resetButton, {
    file_upload$file_upload_sate <- 'reset'
    reset('file')
    })

  FileInput <- reactive({
    if (is.null(file_upload$file_upload_sate)) {
      return(NULL)
    } else if (file_upload$file_upload_sate == 'uploaded') {
      return(input$file)
    } else if (file_upload$file_upload_sate == 'reset') {
      return(NULL)
    }
  })
  
  return(FileInput)
}