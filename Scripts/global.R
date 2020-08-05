#Global File for All_For_One_v1.2

#---Libraries----

library( base )
base_append <- base::append
library( shiny )
library( readr )


#---------------------------------------
# CompletePATH = dirname(sys.frame(1)$ofile)
CompletePATH = paste(getwd(), '/Scripts', sep = '') #-Server
setwd( CompletePATH )

#------Setting up work environment-----
# appPATH = getwd()
# 
# if (substring( appPATH, nchar( appPATH )-8, nchar( appPATH )) != 'CrisPRAL' ){
#  system( 'find /home/ -name "CrisPRAL" > PATH.txt' )
#  
#  file_lines = read_lines( 'PATH.txt' )
#  if (length( file_lines ) >= 2){
#   system( paste ('zenity --warning --width 300 --height 100 
#                   --text="More than one directories named -CrisPRAL- have been 
#                   found, make sure there is only one"' ))
#   stopApp()
#   stop()
#  }
#  
#  appPATH = read_file( 'PATH.txt' )
#  system( 'rm PATH.txt' )
#  
#  CompletePATH = paste(substring( appPATH, 1, nchar( appPATH )-1 ),"/Scripts", sep = '') # nchar - 1 to remove
#  #"/n" due to .txt parsing
#  
#  setwd( CompletePATH )
#  
# }else{ CompletePATH = appPATH }


#-----------Fast loading data-----------
# load(paste(CompletePATH, '/Scripts/Init.Rdata', sep = ''))
#---------------------------------------

#shinyEngineering packages:
shinypckgs <- c("httpuv","shinythemes", "shinyWidgets", "shinydashboard", 'shinyFiles', 'shinyBS', 'shinyjs')

for(package in shinypckgs){
 if (lapply(package, require, character.only = TRUE) == FALSE){
  lapply(package, install.packages, character.only = TRUE)
  lapply(package, require, character.only = TRUE)
 }}

#BiocManager for bioconductor package handling
if (!requireNamespace("BiocManager", quietly=TRUE)){
 install.packages("BiocManager")
 require("BiocManager")
}

#dependant packages required
deppkgs <- c('devtools', 'tidyverse', 'kableExtra', 'rmarkdown', 'processx', 'DT', 'htmlwidgets', 'seqRFLP', 'farver')

for(package in deppkgs){
 
 if (lapply(package, require, character.only = TRUE) == FALSE){
  
  lapply(package, install.packages, character.only = TRUE)
  lapply(package, require, character.only = TRUE)
  
 }
}
count <- dplyr::count
#devtools::install_github() bugfix
options("download.file.method" = "libcurl")

#Bioconductor packages
Biocpkgs <- c("BiocGenerics", "BiocStyle", "BSgenome", 'Biostrings','msa', 'msaR', "checkmate", "ggdendro",
          "reshape2", "Rsamtools", "scales", "ShortRead",  "viridis", "viridisLite", "zoo", "mikelove/fastqcTheoreticalGC", 
          "ngsReports", "wleepang/shiny-directory-input", 'UofABioinformaticsHub/shinyNgsreports')

# BiocManager::install(Biocpkgs)

for (package in Biocpkgs){

 if (lapply(package, require, character.only = TRUE) == FALSE){

  lapply(package, BiocManager::install)
  lapply(package, require, character.only = TRUE)

 }
}

#Github packages
gitpkgs <- c('ropensci/plotly')

for (package in gitpkgs){


 if (lapply(str_split(package, pattern = '/')[[1]][2], require, character.only = TRUE) == FALSE){

  lapply(paste('https://github.com/', package, ".git", sep = '' ), devtools::install_git)
  lapply(str_split(package, pattern = '/')[[1]][2], require, character.only = TRUE)

 }
}

if ( require( 'tinytex' ) == FALSE ){
 
 install.packages( 'tinytex' )
 tinytex::install_tinytex()
 
}

if ( require( 'webshot' ) == FALSE ){
  install.packages( 'webshot' )
  library( webshot )
  webshot::install_phantomjs()
}

#---Modules_Load----

source( str_c( CompletePATH, '/Modules/Functions.R' ))
source( paste( CompletePATH, '/Modules/input-reset_module.R', sep = '' ))
source( paste( CompletePATH, '/Modules/downloadFileModule.R', sep = '' ))

#---Preparation for Lazy Filter----

if ( isFALSE( require( xlsx ))){install.packages( 'xlsx' )
  require( xlsx )}
require( tidyverse )

#----Custom functions----

press_enter <- '
$(function() {
  var $els = $("[data-proxy-click]");
  $.each(
    $els,
    function(idx, el) {
      var $el = $(el);
      var $proxy = $("#" + $el.data("proxyClick"));
      $el.keydown(function (e) {
        if (e.keyCode == 13) {
          $proxy.click();
        }
      });
    }
  );
});
'
