#Global File for All_For_One_v1

#---Libraries----
library( base )
base_append <- base::append
library( shiny )
library( plotly )
library( tidyverse )

if ( require( 'devtools' ) == FALSE ){
  
  install.packages( 'devtools' )
  library( devtools )
  devtools::install_github( 'ropensci/plotly' )
  library( plotly )
  
}

if ( require( 'kableExtra' ) == FALSE ){
  
  install.packages( 'kableExtra' )
  library( kableExtra )
  
}

if ( require( 'rmarkdown' ) == FALSE ){
  
  install.packages( 'rmarkdown' )
  library( rmarkdown )
  
}
if ( require( 'shinyBS' ) == FALSE ){
  
  install.packages( 'shinyBS' )
  library( 'shinyBS' )
  
}

if ( require( 'webshot' ) == FALSE ){
  
  install.packages( 'webshot' )
  library( webshot )
  webshot::install_phantomjs()
  
}

if ( !require( "processx" )) install.packages( "processx" )

#To install Biostrings:
# source( "https://bioconductor.org/biocLite.R" )
# biocLite( "Biostrings" )

library( Biostrings ) 
library( DT )

#------Setting up work environment-----

appPATH = getwd()

if (substring( appPATH, nchar( appPATH )-8, nchar( appPATH )) != 'CrisPRAL' ){
    system( 'find /home/ -name "CrisPRAL" > PATH.txt' )
  
    file_lines = read_lines( 'PATH.txt' )
    if (length( file_lines ) >= 2){
      system( paste ('zenity --warning --width 300 --height 100 
                  --text="More than one directories named -CrisPRAL- have been 
                  found, make sure there is only one"' ))
   stopApp()
   stop()
    }
  
  appPATH = read_file( 'PATH.txt' )
  system( 'rm PATH.txt' )
  
  CompletePATH = substring( appPATH, 1, nchar( appPATH )-1 ) # nchar - 1 to remove 
  #"/n" due to .txt parsing
  
  setwd( CompletePATH )
  
}else{ CompletePATH = appPATH }

#---Modules_Load----

source( paste( CompletePATH, '/Modules/Functions.R', sep = '' ))
source( paste( CompletePATH, '/Modules/input-reset_module.R', sep = '' ))
source( paste( CompletePATH, '/Modules/downloadFileModule.R', sep = '' ))

#---Preparation for Lazy Filter----

if ( isFALSE( require( xlsx ))){install.packages( 'xlsx' )
  require( xlsx )}
require( tidyverse )
