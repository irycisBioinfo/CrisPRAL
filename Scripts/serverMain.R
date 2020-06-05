#All For One main server - Combines all mosaic script analyser and Lazy panel filtering.

serverMain <- function(input, output, session) {
 #Uncomment for app to stop when session if closed:
 
 # session$onSessionEnded(function() {
 #   stopApp()
 # })
 
 options( shiny.maxRequestSize = 100 * 1024 ^ 2 )
 
 #---Mosaic Analyser SERVER----
 
 source( './serverMosaic.R', local = TRUE )
 
 #----Run-Pipeline SERVER----
 
 source( './serverRunPipe.R', local = TRUE )
 
 #----Run-Pipeline sample SERVER----
 
 source( './serverRunExample.R', local = TRUE )
 
 #----Load Alignemnt SERVER----
 
 source( './serverAlignment.R', local = TRUE )
 
 #----Mosaic Grapher SERVER----
 
 source( './serverGraphs.R', local = TRUE )
 
 #---Lazy Panel SERVER----
 
 source( './serverLazy.R', local = TRUE )
 
 #----Fastqc SERVER----
 
 source( './serverFastqc.R', local = TRUE)
 
 
}