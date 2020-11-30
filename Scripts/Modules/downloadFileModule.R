#Generates functions for the automatic creation of an UI and a server for
#downloadable files.

downloadButtonModule <- function(id, label = ''){
  
  ns <- NS(id)
  
  downloadButton(ns('download'), label)
}

downloadCompressedFile <- function(input, output, session, name){
 
 ns <- session$ns
 
 output$download <- downloadHandler(
  filename = function() {
   paste('fastqc_reports','gz', sep = '.')
  },
  content = function(file) {
   
   file <- file.copy(paste(name(),'gz', sep = '.'), file)
   
  }
 )
 
}

downloadCSV <- function(input, output, session, data, name){
  
  ns <- session$ns
  
  output$download <- downloadHandler(
    filename = function() {
      paste(name(),'csv', sep = '.')
    },
    content = function(file) {
      
      data1 <- data() %>% rownames_to_column() %>% rename(Tag = 'rowname')
      file <- write_csv(data1, file)
      
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

# downloadMSA <- function(input, output, session, data, name){
#  
#  ns <- session$ns
#  
#  output$download <- downloadHandler(
#   filename = function(){
#    paste(name(), 'pdf', sep = '.')
#   },
#   content = function(file) {
#    
#    out <- msaPrettyPrint(data(), output = "pdf", showNames="none",
#                          consensusColor="ColdHot", askForOverwrite = FALSE, 
#                          showLogo="top", showLegend = FALSE, 
#                          furtherCode=c("\\textbf{Multiple Sequence Alignment: ClustalW} ",
#                                        "\\DNAgroups{GAR,CTY}","\\defconsensus{.}{lower}{upper}",
#                                        "\\showruler{1}{top}"))
#    file.rename(out, file)
#    
#   })
# }
