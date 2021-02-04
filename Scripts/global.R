#Global File for All_For_One_v1.2

#---Base Libraries----

library( base )
base_append <- base::append
library( shiny )
library( readr )
library("BiocManager")

#---Preparation for Lazy Filter----

require( xlsx )
require( tidyverse )

#---------------------------------------

# CompletePATH = dirname(sys.frame(1)$ofile)
CompletePATH = paste(getwd(), sep = '')
ScriptsPATH = paste(getwd(), '/Scripts', sep = '')
setwd( CompletePATH )

#shinyEngineering packages:
shinypckgs <- c("shinythemes", "shinyWidgets", "shinydashboard", 'shinyFiles', 'shinyBS', 'shinyjs','rintrojs', 'shinycssloaders')
lapply(shinypckgs, require, character.only = TRUE)

#dependant packages required
deppkgs <- c('devtools', 'tidyverse', 'kableExtra', 'rmarkdown', 'processx', 'DT', 'htmlwidgets', 'seqRFLP', 'farver','R.utils')
lapply(deppkgs, require, character.only = TRUE)
# Fixing some masked funcitons:
reset <- shinyjs::reset
count <- dplyr::count
#devtools::install_github() bugfix
options("download.file.method" = "libcurl")

#Bioconductor packages
Biocpkgs <- c("BiocGenerics", "BiocStyle", "BSgenome", 'Biostrings','msa', 'msaR', "checkmate", "ggdendro",
          "reshape2", "Rsamtools", "scales", "ShortRead",  "viridis", "viridisLite", "zoo", "mikelove/fastqcTheoreticalGC", 
          "ngsReports", "wleepang/shiny-directory-input", 'UofABioinformaticsHub/shinyNgsreports')
lapply(Biocpkgs, require, character.only = TRUE)


#Github packages
gitpkgs <- c('ropensci/plotly')

for (package in gitpkgs){

 if (lapply(str_split(package, pattern = '/')[[1]][2], require, character.only = TRUE) == FALSE){

  lapply(paste('https://github.com/', package, ".git", sep = '' ), devtools::install_git)
  lapply(str_split(package, pattern = '/')[[1]][2], require, character.only = TRUE)

 }
}

require( tinytex )
# tinytex::install_tinytex()

# Installation of google-chrome (.deb) is required as well as chromium-browser (apt-get)
require( webshot2 )

#---Modules_Load----

source( str_c( ScriptsPATH, '/Modules/Functions.R' ))
source( str_c( ScriptsPATH, '/Modules/input-reset_module.R'))
source( str_c( ScriptsPATH, '/Modules/downloadFileModule.R'))
source( str_c( ScriptsPATH, '/Modules/are_indels_in_range.R'))

#----Custom functions----

#### ALLOWS USING KEYBORAD ENTER AS A MOUSE CLICK ####.

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
