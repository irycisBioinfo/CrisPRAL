#!/usr/bin/Rscript

filter.table <- function(Tabla, name = 'out' ,path = './',params.mm = 0, params.indel = 0, params.score = 0, params.range = NULL){

#### DESCRIPTION: #####
#
# This script filters csv tables given certain parameters and generates multiple
# csv's with the intermediate application of those filters.
# 
# For the moment parameters are fixed at 0 mismatches, 0 indels.
# In addition if params.score = 'detect', a scoring threshold will be chosen
# based on the scoring distribution. ( This is usefull in the cases where 
# there are no negative scores despite there being some very bad alignments.
# This can happen when aligning using 'overlap' instead of 'global').
#
#######################

library(tidyverse)
library(openxlsx)

# Make sure path ends with a / character.
  
if(str_sub(path, start = nchar(path)) != '/'){
  path <- str_c(path, '/', sep = '')
}
# Mismatch threshold
params.mm <- params.mm
# indel Threshold
params.indel <- params.indel
# score Threshold
if(params.score == 'detect'){
  # We create a leading column and subtract the base data to the lead in order 
  # to obtain the difference between each subsequent value. Max difference obtained
  # and used to infer the threshold.

  mod_table <- Tabla %>% arrange(score) %>% select(score) %>% 
    mutate(scorelaged = lead(score)) %>% mutate(diffScore = scorelaged - score) #%>% 
    #drop_na() %>% summarise(max(diffScore)) %>% pull()
  params.score <- (mod_table$score[which.max(mod_table$diffScore)]+mod_table$score[which.max(mod_table$diffScore)+1])/2
  
}
# score Threshold if not 'detect'
params.score <- params.score

data.filtered.score_threshold <- Tabla %>% filter(score <= params.score)
if(!is.null(params.range)){
  # When writing boolean conditions on an excel, TRUE and FALSE transform to 1 & 0s.
  data.filtered.score_threshold <- data.filtered.score_threshold %>% mutate(overlaps_with_range = case_when(!overlaps_with_range ~ 'no', TRUE ~ 'yes'))
  }
# We need to keep data.filtered.indel though as TRUE & FALSE to filter more efficiently.
data.filtered.score <- Tabla %>% filter(score > params.score) %>% mutate(Freq = 100*Abundance/sum(Abundance))
if(!is.null(params.range)){
  # When writing boolean conditions on an excel, TRUE and FALSE transform to 1 & 0s.
  data.filtered.score <- data.filtered.score %>% mutate(overlaps_with_range = case_when(!overlaps_with_range ~ 'no', TRUE ~ 'yes'))
}
data.filtered.indel <- data.filtered.score %>% filter(mismatch == params.mm) %>% filter(Deletions > params.indel | Insertions > params.indel)

data.filtered.mm <- data.filtered.score %>% filter(mismatch > params.mm) %>% filter(Deletions == params.indel) %>% filter(Insertions == params.indel)
data.filtered.mm_indel <- data.filtered.score %>% filter(mismatch > params.mm) %>% filter(Deletions > params.indel | Insertions > params.indel)
data.filtered.wt <- data.filtered.score %>% filter(mismatch == params.mm) %>% filter(Deletions == params.indel & Insertions == params.indel)


if( !is.null(params.range) ){
  # When writing boolean conditions on an excel, TRUE and FALSE transform to 1 & 0s.
  data.filtered.indel.inrange <- data.filtered.indel %>% filter( overlaps_with_range == 'yes')
  data.filtered.indel.notinrange <- data.filtered.indel %>% filter( overlaps_with_range == 'no' )
  
}

# Lastly modify Tabla
if( !is.null(params.range) ){
  Tabla <- Tabla %>% mutate(overlaps_with_range = case_when(!overlaps_with_range ~ 'no', TRUE ~ 'yes'))
}

# Statistics:

total_reads <- Tabla %>% summarise(total_reads = sum(Abundance))
score <- data.filtered.score_threshold %>% summarise(threshold_score = sum(Abundance))

wt <- data.filtered.wt %>% summarise(wildtypes = sum(Abundance))
mm <- data.filtered.mm %>% summarise(mismatches = sum(Abundance))
indels <- data.filtered.mm_indel %>% summarise("indels with mismatches" = sum(Abundance))


if( !is.null(params.range) ){
  clean_indels_inrange <- data.filtered.indel.inrange %>% summarise("indels in range" = sum(Abundance))
  clean_indels_notinrange <- data.filtered.indel.notinrange %>% summarise("indels not in range" = sum(Abundance))
}else{
  clean_indels <- data.filtered.indel %>% summarise(indels = sum(Abundance))
}


if( !is.null(params.range) ){
  statistics <- bind_cols(total_reads, score, wt,mm,indels, clean_indels_notinrange, clean_indels_inrange)
}else{
  statistics <- bind_cols(total_reads, score, wt,mm,indels, clean_indels)
}
# Creating Excel
wb <- createWorkbook()

openxlsx::addWorksheet(wb, 'MAIN')
writeDataTable(wb, Tabla, sheet = 'MAIN')

openxlsx::addWorksheet(wb, paste('threshold_score -', params.score , sep = ' '))
writeDataTable(wb, data.filtered.score_threshold, sheet = paste('threshold_score -', params.score , sep = ' '))

openxlsx::addWorksheet(wb, 'wildtype')
writeDataTable(wb, data.filtered.wt, sheet = 'wildtype')

openxlsx::addWorksheet(wb, 'mismatches')
writeDataTable(wb, data.filtered.mm, sheet = 'mismatches')

openxlsx::addWorksheet(wb, 'indels_with_mismatches')
writeDataTable(wb, data.filtered.mm_indel, sheet = 'indels_with_mismatches')

if(!is.null(params.range)){
  
  openxlsx::addWorksheet(wb, paste0('indels_not_in_range - ',params.range[1],'-',params.range[length(params.range)]))
  writeDataTable(wb, data.filtered.indel.notinrange, sheet = paste0('indels_not_in_range - ',params.range[1],'-',params.range[length(params.range)]))
  
  openxlsx::addWorksheet(wb, paste0('indels_in_range - ',params.range[1],'-',params.range[length(params.range)]))
  writeDataTable(wb, data.filtered.indel.inrange, sheet = paste0('indels_in_range - ',params.range[1],'-',params.range[length(params.range)]))

}else{
  openxlsx::addWorksheet(wb, 'clean_indels')
  writeDataTable(wb, data.filtered.indel, sheet = 'clean_indels')
}

openxlsx::addWorksheet(wb, 'statistics')
writeDataTable(wb, statistics, sheet = 'statistics')

return(c(workbook = wb))

}
