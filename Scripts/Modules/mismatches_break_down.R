
library(tidyverse)
library(Biostrings)

mismatches_break_down <- function(alignment_raw, data){
  
  #### DESCRIPTION: #####
  #
  # This script breaks down the position and the change of all mismatches 
  # found in the data.
  #
  # alignment_raw is the results of applying Clusters_Alignments() from Functions.R
  # 
  #
  #######################

mmTable <- mismatchTable(alignment_raw$aln)
mmTable <- mmTable %>% select(-PatternStart, -PatternEnd, SubjectEnd)
mmTable <- mmTable %>% group_by(PatternId)
mmTable <- mmTable %>% group_by(PatternId) %>% mutate(SNP = str_c(SubjectStart,SubjectSubstring,'>',PatternSubstring))
mmTable <- mmTable %>% select(1,6) %>% mutate(SNP = str_c(SNP, collapse = ';')) %>% ungroup() %>% distinct()

Unsort_Tabla_SNP <- left_join(alignment_raw$tmp %>% rownames_to_column() %>% 
                                mutate('rowname' = as.integer(rowname)), 
                                mmTable %>% rename(rowname = PatternId), by = 'rowname') %>% 
                                column_to_rownames('rowname') %>% select(ID,SNP)
final_table <- left_join(data,Unsort_Tabla_SNP, by = 'ID')

return(final_table)

}