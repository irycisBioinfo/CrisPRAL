#!/usr/bin/Rscript

################################################################################
############################ Isoforma corta fixer ##############################
#                                                                              #
# Description:                                                                 #
#                                                                              #
# From a VCF file, substracts a fixed amount to the position column           #
#                                                                              #
# Usage:                                                                       #
#                                                                              #
# fix_pos_isoforma_corta.R *.vcf [INT]                                         #
#                                                                              #
################################################################################

##################
# TODO: asserts  #
##################

library(tidyverse)

f <- commandArgs(trailingOnly=TRUE)
# substract 118 to every line in POS column of vcf
correction.int <- as.numeric(grep(f, invert = TRUE, value = TRUE, pattern = '\\.'))

for(file in grep(f, pattern = '.vcf', value = TRUE)){
  
  vcf <- read_tsv(file, comment = "##", col_names = TRUE)
  vcf <- vcf %>% mutate(POS = as.numeric(POS) + correction.int)
  
  system(paste('cat', file, '| grep ^## >', '`echo', file, '| sed s/.vcf/.fix.vcf/`', sep = ' '))
  
  write_tsv(vcf, str_replace(file, pattern = ".vcf", replacement = ".fix.vcf"), append = TRUE, col_names = TRUE)
  # Possibly fixed?
  # print("Commas issue not fixed, execute sed s/^\"//g | sed s/\"$//g in your command line to fix it")
}