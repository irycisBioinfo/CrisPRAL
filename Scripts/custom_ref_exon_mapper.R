#!/usr/bin/Rscript

################################################################################
####################### Custom reference exon mapper ###########################
#                                                                              #
# Description:                                                                 #
#                                                                              #
# From a VCF file, maps the position to an exon given in a custom exon -       #
# position table.                                                              #
#                                                                              #
# Usage:                                                                       #
#                                                                              #
# custom_ref_exon_mapper.R [VCF FILE/S] [ANNOT FILE]                           #
#                                                                              #
# Annot file must have a .tsv extension.                                       #  
#                                                                              #
################################################################################

##################
# TODO: asserts  #
##################


library(tidyverse)

f <- commandArgs(trailingOnly=TRUE)

gene_annotation <- read_tsv(grep(f, pattern = '.tsv', value = TRUE), col_names = TRUE)

for(file in grep(f, pattern = '.vcf', value = TRUE)){

vcf <- read_tsv(file, comment = "##", col_names = TRUE)
vcf_pos <- vcf$POS

vcf_chrom <- c()
for(entry in vcf_pos){
  
  exon.idx <- last(which(entry > gene_annotation$INICIO))
  exon <- gene_annotation$EXON[exon.idx]
  vcf_chrom <- c(vcf_chrom, exon)
  
}

new_vcf <- bind_cols(`#CHROM` = vcf_chrom, vcf %>% select(-`#CHROM`))

system(paste('cat', file, '| grep ^## >', '`echo', file, '| sed s/.vcf/.exons.vcf/`', sep = ' '))

write_tsv(new_vcf, str_replace(file, pattern = ".vcf", replacement = ".exons.vcf"), 
          append = TRUE, col_names = TRUE)
}

