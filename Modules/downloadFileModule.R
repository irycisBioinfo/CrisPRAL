#Generates functions for the automatic creation of an UI and a server for
#downloadable files.

downloadButtonModule <- function(id, label = ''){
  
  
  ns <- NS(id)
  
  downloadButton(ns('download'), label)
}


downloadCSV <- function(input, output, session, data, name){
  
  ns <- session$ns
  
  output$download <- downloadHandler(
    filename = function() {
      paste(name(),'csv', sep = '.')
    },
    content = function(file) {
      
      file <- write_csv(data, file)
      
    }
  )}

downloadFASTA <- function(input, output, session, data, name, 
                          indices = FALSE){
  
  ns <- session$ns
  
  if ( isFALSE( indices )){
    output$download <- downloadHandler(
      filename = function(){
        paste(name(), 'fasta', sep = '.')
        },
      content = function(file) {
        
        file <- writeXStringSet(data(), file)
      })
    
  }else{
    output$download <- downloadHandler(
      filename = function(){
        paste(name(), 'fasta', sep = '.')
      },
      content = function(file) {
        
        file <- writeXStringSet(data()[indices()], file)
      })
  }
}

downloadPAIR <- function(input, output, session, data, name){
  
  ns <- session$ns
  
  output$download <- downloadHandler(
    filename = function(){
      paste(name(), 'pair', sep = '.')
    },
    content = function(file) {
      
      file <- writePairwiseAlignments(data(), file)
      
    })
}
