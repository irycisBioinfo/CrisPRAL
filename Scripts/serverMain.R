#All For One main server

############################### DESCRIPTION ##############################
#Combines all mosaic script analyser and Lazy panel filtering.
##########################################################################
serverMain <- function(input, output, session) {
 #Uncomment for app to stop when session if closed:
 
 # session$onSessionEnded(function() {
 #   stopApp()
 # })
 
 # REACTIVE VALUES ACCESIBLE BY THE WHOLE APP #
 datos <- reactiveValues()
 datos$system_administrator <- "sergio.fern1994@gmail.com"
        
 options( shiny.maxRequestSize = 500 * 1024 ^ 2 )
 
 #---Mosaic Analyser SERVER----
 
 source( 'Scripts/serverMosaic.R', local = TRUE )
        
#---BATCH Mosaic Analyser SERVER----

source( 'Scripts/serverMosaic_batch.R', local = TRUE )
 
 #----Run-Pipeline SERVER----
 
 source( 'Scripts/serverRunPipe.R', local = TRUE )
        
 #---- Short Variant Caller SERVER----
        
 source( 'Scripts/serverShortVariantCaller.R', local = TRUE )
 
 #----Run-Pipeline sample SERVER----
 
 source( 'Scripts/serverRunExample.R', local = TRUE )
 
 #----Load Alignemnt SERVER----
 
 source( 'Scripts/serverAlignment.R', local = TRUE )
 
 #----Mosaic Grapher SERVER----
 
 source( 'Scripts/serverGraphs.R', local = TRUE )
 
 #---Lazy Panel SERVER----
 
 source( 'Scripts/serverLazy.R', local = TRUE )
 
 #--- Simple vcf viwer SERVER----
        
 source( 'Scripts/serverSimple_vcf_viewer.R', local = TRUE)
        
observeEvent(input$btn_landing,{return(browser())})
 
 
}