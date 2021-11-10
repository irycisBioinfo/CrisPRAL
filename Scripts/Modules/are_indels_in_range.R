library(tidyverse)
library(vctrs)
library(xlsx)

are_indels_in_range <- function(data, position_range){
  
  #### DESCRIPTION: #####
  #
  # This script finds whether the data table provided has indels falling 
  # within a given range.
  #
  #######################

#######################.
## Set parameters to work with (position range) ----
#######################.

#######################.
## Operations ----
#######################.

##########.
# FIRST: I need to know the maximun number of indel there are.
##########.
max_ndel = data %>% select(contains('deletion_range')) %>% mutate(ndel = str_count(deletion_range, ';')+1) %>% drop_na() %>% summarise(max(ndel)) %>% pull()
max_nin = data %>% select(contains('insert_range')) %>% mutate(nin = str_count(insert_range, ';')+1) %>% drop_na() %>% summarise(max(nin)) %>% pull()

##########.
# SECOND: we need the new column names for the indels
##########.
del_columns = str_c('Del',c(1:max_ndel), sep = '')
in_columns = str_c('Ins',c(1:max_nin), sep = '')

dels_desglosados <- data %>% select(contains('deletion_range')) %>% separate(col = 'deletion_range',into = del_columns,extra = 'merge',fill = 'right', sep = ';') %>%
  mutate( across( everything(), ~replace_na(.,'')) )
ins_desglosados <- data %>% select(contains('insert_range')) %>% separate(col = 'insert_range',into = in_columns,extra = 'merge',fill = 'right', sep = ';') %>% 
  mutate( across( everything(), ~replace_na(.,'')) )

##########.
# THIRD: separate() columns again but this time by ':', into colnames(insert1_from, insert1_to, isnert2_from....)
##########.

into <- str_c(vec_rep_each(in_columns, 2),c("_from","_to")) # Vector prepared for "into" option in separate()

# Inserts:
ins = c()
for(i in c(1:ncol(ins_desglosados))){
  
  into <- str_c(vec_rep_each(in_columns[i], 2),c("_from","_to")) # Vector prepared for "into" option in separate()
  tmp <- ins_desglosados[i] %>% drop_na() %>% separate(col = names(ins_desglosados[i]), into = into, extra = "merge", fill = 'right', sep = ':')
  ins = c(ins, tmp)
}
final_ins_data <- as_tibble(ins) %>% mutate( across( everything(), ~replace_na(.,'')) )

# Deletions:
dels = c()
for(i in c(1:ncol(dels_desglosados))){
  
  into <- str_c(vec_rep_each(del_columns[i], 2),c("_from","_to")) # Vector prepared for "into" option in separate()
  tmp <- dels_desglosados[i] %>% drop_na() %>% separate(col = names(dels_desglosados[i]), into = into, extra = "merge", fill = 'right', sep = ':')
  dels = c(dels, tmp)
}
final_dels_data <- as_tibble(dels) %>% mutate( across( everything(), ~replace_na(.,'')) )

##########.
# FOURTH: find overlaps with provided interest range.
##########.

# We will obtain an extra column for every original column, with boolean values indicating if it is the range of the values given:
is_ins_data_in_range <- final_ins_data %>% mutate(in_range = across(everything(),~case_when(.%in% position_range ~ TRUE, TRUE ~ FALSE)))
is_del_data_in_range <- final_dels_data %>% mutate(in_range = across(everything(),~case_when(.%in% position_range ~ TRUE, TRUE ~ FALSE)))

# Sum up all the extra columns into a single boolean vector.
are_ins_in_range_summary <- rowSums(is_ins_data_in_range$in_range) != 0
are_dels_range_summary <- rowSums(is_del_data_in_range$in_range) != 0

# Fuse both
indel_in_range_summary <- rowSums(cbind(are_dels_range_summary, are_ins_in_range_summary)) != 0

# Since the order of the data was never altered we can simply bind the new column
data_with_range_check <- cbind(data, overlaps_with_range = indel_in_range_summary)

return(data_with_range_check)
}