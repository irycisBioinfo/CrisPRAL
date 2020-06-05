downloadButtonModule <- function(id, label = ''){
  
  ns <- NS(id)
  
    downloadButton(ns('download'), label)
  }


downloadCSV <- function(input, output, session, data, name){
  
  ns <- session$ns

  output$download <- downloadHandler(
    filename = function() {
      paste(name,'csv', sep = '.')
    },
    content = function(file) {
    
    file <- write_csv(data, file)
  
    }
    )}
