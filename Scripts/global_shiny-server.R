#Global File for All_For_One_server-v0.1.2

#---Libraries----
library( base )
base_append <- base::append
library( shiny )
count <- dplyr::count

#shinyEngineering packages:
shinypckgs <- c("shinythemes", "shinyWidgets", "shinydashboard", 'shinyFiles', 'shinyBS', 'shinyjs')

for(package in shinypckgs){
 
 if (lapply(package, require, character.only = TRUE) == FALSE){
  
  lapply(package, install.packages, character.only = TRUE)
  lapply(package, require, character.only = TRUE)
  
 }}
#BiocManager for bioconducter package handling
if (!requireNamespace("BiocManager", quietly=TRUE)){
 install.packages("BiocManager")
 require("BiocManager")
}

#dependant packages required
deppkgs <- c('devtools', 'tidyverse', 'kableExtra', 'rmarkdown', 'processx', 'DT')

for(package in deppkgs){
 
 if (lapply(package, require, character.only = TRUE) == FALSE){
  
  lapply(package, install.packages, character.only = TRUE)
  lapply(package, require, character.only = TRUE)
  
 }}

count <- dplyr::count
#devtools::install_github() bugfix
options("download.file.method" = "libcurl")

#Bioconductor packages
Biocpkgs <- c("BiocGenerics", "BiocStyle", "BSgenome", 'Biostrings','msa', 'msaR', "checkmate", "ggdendro",
              "reshape2", "Rsamtools", "scales", "ShortRead",  "viridis", "viridisLite", "zoo")

for (package in Biocpkgs){
 
 if (lapply(package, require, character.only = TRUE) == FALSE){
  
  lapply(package, BiocManager::install)
  lapply(package, require, character.only = TRUE)
  
 }}

#packages for fastqc analysis module
gitpkgs <- c('ropensci/plotly', "mikelove/fastqcTheoreticalGC", "UofABioinformaticsHub/ngsReports", 'UofABioinformaticsHub/fastqcRShiny')

for (package in gitpkgs){
 
 if (lapply(str_split(package, pattern = '/')[[1]][2], require, character.only = TRUE) == FALSE){
  
  lapply(paste('https://github.com/', package, ".git", sep = '' ), devtools::install_git)
  # lapply(gitpkgs, devtools::install_github('mikelove/fastqcTheoreticalGC'))
  lapply(str_split(package, pattern = '/')[[1]][2], require, character.only = TRUE)
  
 }}

if ( require( 'tinytex' ) == FALSE ){
 
 install.packages( 'tinytex' )
 tinytex::install_tinytex()
 
}

if ( require( 'webshot' ) == FALSE ){
 install.packages( 'webshot' )
 library( webshot )
 webshot::install_phantomjs()
}

#------Setting up work environment-----

appPATH = getwd()
CompletePATH <- '/home/sergio/ShinyApps/AllForOne/CrisPRAL'

if (appPATH!=CompletePATH){ 
  
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
