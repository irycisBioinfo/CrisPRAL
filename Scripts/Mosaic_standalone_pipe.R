#!/usr/bin/Rscript

#### DESCRIPTION: #######################################################
#
# This scripts is a standalone version of Mosaic's shiny app, allows for cluster
# generation from a sample provided.
#
# Includes some additional features like identification and printing of the
# positions of indels in an auxiliary column.
#
# This script was thought to be used in par with "filter_table.R" 
# and "are_indels_in_Range.R" in scripts/
# 
#### USAGE: ####
#
# ./Mosaic_standalone_pipe.R [R1_FILES] -reference [FILE] -adapters [FILE] -primers [FILE] [OPTIONS] 
#
# OPTIONS:
#
# -adapters [FILE]
# -primers [FILE]
# 
# -output [PATH]                              Location where results will be saved. (default: PATH = './')
# -score_threshold <[INT] | 'detect'>         Final output will have a table with all data with a score below the threshold equal to INT (default: INT = 0)
#                                             If -score_threshold = 'detect', a score is selected based on the scoring distribution.
# -mismatches                                 Flag if present mismatches will be broken down and printed in another column together with the rest of results
# -indel_interest_range <[INT-INT] | NULL>    Range of interest for indels identification. (default: NULL) 
# 
#########################################################################.
################################ SETUP ##################################.
#########################################################################.

library(tidyverse)
library(Biostrings)
library(vctrs)
library(xlsx)

CompletePATH = paste(getwd(), sep = '') #-Server

source(str_c(CompletePATH,'/Modules/aid_functions.R', sep = ''))
source(str_c(CompletePATH,'/Modules/Functions.R', sep = ''))
source(str_c(CompletePATH,'/Modules/filter_table.R', sep = ''))
source(str_c(CompletePATH,'/Modules/are_indels_in_range.R', sep = ''))
source(str_c(CompletePATH,'/Modules/mismatches_break_down.R', sep = ''))

rename <- dplyr::rename
filter <- dplyr::filter

f <- commandArgs(trailingOnly = TRUE)
# f <- c('/media/bioinfo/Seagate_Expansion_Drive/ALMUDENA_JUL_2020/C_S35_L001_R1_001.fastq', ',-reference', 'TYR\ reference.fasta', '-primers', 'tyr_primers.fasta', '-output', 'results_Almudena2/', '-score_threshold',"detect", '-indel_interest_range', '70-170')
files <- grep(f, pattern = 'fastq$', value = TRUE)

# Find primers
primers_pos <- grep(f, pattern = '-primers')
if(!is_empty(primers_pos)){
  primers <- f[primers_pos + 1]
  # Read primers:
  store.primers <- read.primers(primers)
}else{store.primers <- c(F_Primer1 = 'Empty', F_Primer2 = 'Empty', R_Primer1 = 'Empty', R_Primer2 = 'Empty')}

# Find Adapters
adapters_pos <- grep(f, pattern = '-adapters')
if(!is_empty(adapters_pos)){
  adapters <- f[adapters_pos + 1]
  # Read adapters:
  if( is.null(adapters) ){store.adapters <- c(F_adapter1 = 'Empty', R_adapter1 = 'Empty')}
  adapters.lines <- read_lines(adapters)
  store.adapters <- c(F_adapter1 = adapters.lines[2], R_adapter1 = adapters.lines[5])
}else{store.adapters <- c(F_adapter1 = 'Empty', R_adapter1 = 'Empty')}

# Find reference
reference_pos <- grep(f, pattern = '-reference')
if(!is_empty(reference_pos)){
  reference <- f[reference_pos + 1]
}else{print(errorCondition(message = "mandatory -reference option not found"))}

# Find output dir
output_pos <- grep(f, pattern = '-output')
if(!is_empty(output_pos)){
output <- f[output_pos + 1]
}else{output <- './'}

# Find threshold score
score_threshold_pos <- grep(f, pattern = '-score_threshold')
if(!is_empty(score_threshold_pos)){
  score_threshold <- f[score_threshold_pos + 1]
}else{
  score_threshold <- 0
}

# check mismatches flag
mismatches <- grep(f, pattern = '-mismatches')
if(!is_empty(mismatches)){mismatches <- TRUE}else{mismatches <- FALSE}

# Find indel interest range
interest_range_pos <- grep(f, pattern = '-indel_interest_range')
if(!is_empty(interest_range_pos)){
  
  interest_range <- f[interest_range_pos + 1]

  interest_range <- as.numeric(str_split(interest_range, pattern = '-')[[1]])
  interest_range <- c(interest_range[1]:interest_range[2])
  
}else{
  interest_range <- NULL
}

# Make sure output ends with a / character.
if(str_sub(output, start = nchar(output)) != '/'){
  output <- str_c(output, '/', sep = '')
}

# Directory handling:

if(!dir.exists(output)){dir.create(output)}

for(file_pair1 in files){
  
    file_id <- str_split(file_pair1, pattern = '/')[[1]][length(str_split(file_pair1, pattern = '/')[[1]])]
    name <- str_split(file_id, pattern = '_R1_')[[1]][1]
    tmppipelinedir <- str_c(output, name, sep = '')
  
    if(!dir.exists(tmppipelinedir)){dir.create(tmppipelinedir)}
    
    # Reverse file identification:
    
    file_pair2 <- str_replace_last(file_pair1, pattern = 'R1', replacement = 'R2')
    
  ###########################.
  ## EXECUTE ----
  ###########################.
  
    command = paste0(
      CompletePATH,
      "/pipeline_vUMI.pl --r1 ",
      file_pair1,
      " --r2 ",
      file_pair2,
      " --single_end ",
      'FALSE',
      " --min_len ",
      100,
      " --Adapter_R1 ",
      store.adapters['F_adapter1'],
      " --Adapter_R2 ",
      store.adapters['R_adapter1'],
      " --trimA1 ",
      0,
      " --trimA2 ",
      0,
      " --F_Primer1 ",
      store.primers['F_Primer1'],
      " --F_Primer2 ",
      store.primers['F_Primer2'],
      " --R_Primer1 ",
      store.primers['R_Primer1'],
      " --R_Primer2 ",
      store.primers['R_Primer2'],
      " --trimP1 ",
      0,
      " --trimP2 ",
      0,
      " --cov ",
      1,
      " --id ",
      1,
      " --tmpdir ",
      as.character(tmppipelinedir),
      collapse = ""
    )
    
    system(command)
  
  if(!file.exists(paste(tmppipelinedir,"/cluster.tsv", sep = ''))){
    #-assert no missing files
    
    print('::ERROR:: Pipeline was broken; check files, extensions, reload and if problem persists contact administrator at sergio.fern1994@gmail.com')
    
  }
  
  #This directory has to be reactive
  fasta <- read_tsv(paste(tmppipelinedir,"/cluster.tsv", sep = ''), 
                          col_names = FALSE)
  cluster <- read_tsv(paste(tmppipelinedir,"/cluster.bak.tsv", sep = ''), 
                            col_names = FALSE)
  
  if (isEmpty(colnames(fasta))){#-If dataset is completely filtered 
  
    print("::WARNING:: all data was filtered.")
  }else{
    colnames(fasta) <- c("ID", "Sequence")
    colnames(cluster) <- c("ClusterN", "Length", "ID", "Identity")
  } 
  
  allfasta <- read_delim(paste(tmppipelinedir,"/final.fasta" ,sep = ''), 
                               delim = ' ', col_names =FALSE) %>% select(X1)
  
  justfasta <- allfasta %>% filter(!str_starts(X1, '>'))
  justIds <- allfasta %>% filter(str_starts(X1, '>'))
  allfasta <- bind_cols(justIds, justfasta)
  colnames(allfasta) <- c('ID','SEQ')
  colnames(cluster) <- c("ClusterN", "Length", "ID", "Identity")
  #remove '>' character of fasta identifier
  allfasta <- allfasta %>% mutate_at(vars(ID), list(~if_else(str_starts(., '>'),
                                                              true = str_extract(., '(?<=>)[:graph:]+'), 
                                                              false = .)))
  
  clusterF <- cluster %>% select("ClusterN", "ID", "Identity")
  colnames(clusterF) <- c('CLUSTER', 'ID', 'Identity')
  clufast <- full_join(allfasta, clusterF) %>% group_by(CLUSTER) %>% nest() %>% arrange(CLUSTER)
  
  #Tabla_raw is a rearrangement of cluster to prepare the data for
  #fusing with Cluster-Sequences aswell as integrating with Alignment data
  
  Tabla_raw = cluster %>%
    group_by(ClusterN) %>%
    mutate(Abundance = n()*2) %>%
    ungroup() %>%
    mutate(Freq = 100 *Abundance / (n()*2)) %>%
    filter(Identity == "*") %>%
    select(ID, Abundance, Freq) %>%
    arrange(desc(Abundance))
  
  #----Associating clusterN with sequences----
  
  Tabla_raw_clstrs <- Tabla_raw %>% rowid_to_column('ClusterN') %>% 
    select(-Abundance, -Freq)
  
  clustering <- full_join(Tabla_raw_clstrs, clufast %>% 
                                  unnest(), by = 'ID') %>% arrange(CLUSTER)
  clustering <- clustering %>% fill(ClusterN, .direction = 'down') %>% 
    select(-Identity, -CLUSTER ) %>% 
    group_by(ClusterN) %>% nest %>% 
    arrange(ClusterN)
  
  #Preparing for Alignment operations---- 
  #-asserting Reference does not provide an error:
  
  Sequences <- readDNAStringSet(paste(tmppipelinedir,"/cluster", sep = ''))
  Ref <- readDNAStringSet(reference)
  #Alignment for Sequences with Reference
  ls_ClusterAln_Ref = Clusters_Alignments(Sequences, Ref, type = 'overlap') #Function in Functions.R module.
  #Another Tabla variable is declared as an unsorted version together with
  #an associated Abundance for INDEL processing in graphing section.
  Tabla_unsort <- ls_ClusterAln_Ref$tmp %>% rename(Deleted_bps = deletions) %>% rename(Inserted_bps = insertions)
  
  
  aln = ls_ClusterAln_Ref$aln
  
  #Changing deletions and insertion columns to show total number of indels 
  #instead of total number of indel bps
  number_of_dels <- map(as.list(deletion(aln)), ~length(.))
  total_deletions_per_cluster <- as_tibble(purrr::flatten_dbl(number_of_dels)) %>% rename(Deletions = value)
  number_of_ins <- map(as.list(insertion(aln)), ~length(.))
  total_insertions_per_cluster <- as_tibble(purrr::flatten_dbl(number_of_ins)) %>% rename(Insertions = value)
  
  # Tabla_unsort_total_indels <- cbind(cbind(Tabla_unsort %>% select(-Deleted_bps, -Inserted_bps), total_deletions_per_cluster, total_insertions_per_cluster), Tabla_unsort %>% select(Deleted_bps, Inserted_bps))
  
  Tabla_unsort_total_indels <- cbind(Tabla_unsort %>% select(-Deleted_bps, -Inserted_bps), total_deletions_per_cluster, Tabla_unsort %>% select(Deleted_bps), total_insertions_per_cluster, Tabla_unsort %>% select(Inserted_bps))
  
  Tabla = inner_join(Tabla_raw, Tabla_unsort_total_indels) %>% 
    mutate(score = round(score,1), Freq = signif(Freq,2)) %>% 
    arrange(desc(Abundance)) %>% select(-width, -start, -end)
  
  rownames(Tabla) <- str_c('Cluster', rownames(Tabla))
  
  id <- str_split(names(unaligned(aln@pattern)), pattern = ' ')
  
  aligned_names <- c()
  for(i in id){aligned_names <- c(aligned_names, i[1])}
  
  pos_in_table <- match(aligned_names, Tabla$ID)
  
  insertions <- as.data.frame(indel(aln)@insertion) %>% select(-group_name)
  deletions <- as.data.frame(indel(aln)@deletion) %>% select(-group_name)
  
  insertions_pos <- insertions %>% mutate(middle_pos = paste(start,":",end, sep = '')) %>% select(-start,-end,-width)
  deletions_pos <- deletions %>% mutate(middle_pos = paste(start,":",end, sep = '')) %>% select(-start,-end,-width)
  
  insertions_clst <- insertions_pos %>% mutate(cluster = str_c('Cluster',pos_in_table[insertions_pos$group])) %>% 
    select(-group) %>% group_by(cluster) %>% summarise(insert_range = paste(middle_pos, collapse = ";"))
  deletions_clst <- deletions_pos %>% mutate(cluster = str_c('Cluster',pos_in_table[deletions_pos$group])) %>% 
    select(-group) %>% group_by(cluster) %>% summarise(deletion_range = paste(middle_pos, collapse = ";"))
  
  Tabla_with_rownames <- Tabla %>% rownames_to_column('cluster')
  
  Tabla_with_ins <- left_join(Tabla_with_rownames, insertions_clst, by = 'cluster')
  Tabla_with_indels <- left_join(Tabla_with_ins, deletions_clst, by = 'cluster')
  Tabla_final <- Tabla_with_indels
  #################################################.
  ## Identify if indels falls in interest range
  #################################################.
  
  if(!is.null(interest_range)){
    Tabla_final <- are_indels_in_range(Tabla_final, interest_range)
  }
  
  #################################################.
  ## Identify mismatch positioning and change
  #################################################.
  if(mismatches){
    Tabla_final <- mismatches_break_down(ls_ClusterAln_Ref, Tabla_final)
  }
  
  ######################################################################.
  ## Filter final table and structure data into different xlsx sheets
  ######################################################################.
  
  filter.table_results <- filter.table(Tabla = Tabla_final, path = output, name = name, params.score = score_threshold, params.range = interest_range)
  
  # Save data ready for delivery
  openxlsx::saveWorkbook(filter.table_results$workbook, paste0(output, name, '.xlsx'), overwrite = TRUE)
  
}
