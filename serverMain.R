#All For One main server - Combines mosaic analyser and Lazy panel filtering.

# Big Issue: webshot does not work in Linux Debian, must comment out all graph
# printing for pdf output.

serverMain <- function(input, output, session) {
  #Uncomment for app to stop when session is closed:
  
  # session$onSessionEnded(function() {
  #   stopApp()
  # })
  
  options( shiny.maxRequestSize = 100 * 1024 ^ 2 )
  
  #---Mosaic Analyser SERVER----
  
  source( 'serverMosaic.R', local = TRUE )
  
  #---Lazy Panel SERVER----
  
  source( 'serverLazy.R', local = TRUE )

  
}