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
# -output [PATH]                              Location where results will be saved. (default: './')
# -alignment <"global"|"overlap">             Selects alignment method to perform. (default: 'overlap')
# -score_threshold <[INT] | 'detect'>         Final output will have a table with all data with a score below the threshold equal to INT (default: INT = 0)
#                                             If -score_threshold = 'detect', a score is selected based on the scoring distribution.
# -mismatches                                 Flag if present mismatches will be broken down and printed in another column together with the rest of results
# -indel_interest_range <[INT-INT] | NULL>    Range of interest for indels identification. (default: NULL)
# -append_sequences                           Appends representative sequences of the clusters to the table as a final column
# 
#########################################################################.
################################ SETUP ##################################.
#########################################################################.

library(tidyverse)
library(Biostrings)
library(vctrs)
library(xlsx)

CompletePATH = paste(getwd(), sep = '') #-Server

source(str_c(CompletePATH,'/Scripts/Modules/Functions.R', sep = ''))
source(str_c(CompletePATH,'/Scripts/Modules/filter_table.R', sep = ''))
source(str_c(CompletePATH,'/Scripts/Modules/are_indels_in_range.R', sep = ''))
source(str_c(CompletePATH,'/Scripts/Modules/aid_functions.R', sep = ''))

rename <- dplyr::rename
filter <- dplyr::filter

f <- commandArgs(trailingOnly = TRUE)
# f <- c('/media/bioinfo/Seagate_Expansion_Drive/ALMUDENA_JUL_2020/C_S35_L001_R1_001.fastq', ',-reference', 'TYR\ reference.fasta', '-primers', 'tyr_primers.fasta', '-output', 'results_Almudena2/', '-score_threshold',"detect", '-indel_interest_range', '70-170')
files <- grep(f, pattern = 'fastq$', value = TRUE)

######################################################.
###################### OPTIONS #######################
######################################################.

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

# check append_sequences flag
sequences <- grep(f, pattern = '-append_sequences')
if(!is_empty(sequences)){sequences <- TRUE}else{sequences <- FALSE}

# check alignment options
alignment_pos <- grep(f, pattern = '-alignment')
if(!is_empty(alignment_pos)){alignment <- f[alignment_pos + 1]}else{alignment <- 'overlap'}

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

#######################################################.
#######################################################.

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
      CompletePATH,'/Scripts',
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
  
  Tabla_raw = cluster %>%
    group_by(ClusterN) %>%
    mutate(Abundance = n()*2) %>%
    ungroup() %>%
    mutate(Freq = 100 *Abundance / (n()*2)) %>%
    filter(Identity == "*") %>%
    select(ID, Abundance, Freq, Length) %>%
    arrange(desc(Abundance))
  
  #Preparing for Alignment operations---- 
  #-asserting Reference does not provide an error:
  
  Sequences <- readDNAStringSet(paste(tmppipelinedir,"/cluster", sep = ''))
  Ref <- readDNAStringSet(reference)
  #Alignment for Sequences with Reference
  ls_ClusterAln_Ref = Clusters_Alignments(Sequences, Ref, type = alignment) #Function in Functions.R module.
  aln = ls_ClusterAln_Ref$aln

  Tabla_final = inner_join(Tabla_raw, ls_ClusterAln_Ref$tmp) %>% 
    mutate(score = round(score,1), Freq = signif(Freq,2)) %>% 
    arrange(desc(Abundance))# %>% select(-width, -start, -end)
  
  rownames(Tabla_final) <- str_c('Cluster', rownames(Tabla_final))
  
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
  #################################################.
  ## Append sequences as the last column
  #################################################.
  if(sequences){
    # This operations have the objective of being able to compare sequences of different samples between each other.
    Tabla_final <- inner_join(Tabla_final, allfasta)
  }
  
  ######################################################################.
  ## Filter final table and structure data into different xlsx sheets
  ######################################################################.
  
  filter.table_results <- filter.table(Tabla = Tabla_final, path = output, name = name, params.score = score_threshold, params.range = interest_range)
  
  # Save data ready for delivery
  openxlsx::saveWorkbook(filter.table_results$workbook, paste0(output, name, '.xlsx'), overwrite = TRUE)
  
}
