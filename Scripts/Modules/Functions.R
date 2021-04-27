###################.
# DESCRIPTION:
#
# This script stores different multiporpouse functions required
# by different scripts. A function is declared here if the extension of it
# is larger than 5 code lines or if it serves multiple scripts.
#
# Functions:
#
# Position_vector() - used in serverAlignment.R
# 
# lengthFixing() - used in Position_vector()
#
# InDel_Extraction() - used in serverGraphs.R
#
# InDel_Sizes() - used in InDel_Extraction()
#
# Unzip_file() - used in multiple scripts
#
# Unravel_Positions() - used in serverGraphs.R
#
# find_target_location() - used in serverMosaic.R
#
# Group_list() - used in serverGraphs.R
#
# Clusters_Alignments() - used in serverRunPipe.R
#
# file.dir() - used in serverMosaic.R and others
#
# conditional_input_decision_tree() - used in serverMosaic_batch.R
###################.

Start_gaps_size <- function(string){
  
  ##################
  # Description:
  #
  # Returns the size of the gaps introduced in the start of an alignment
  # represented by dashes.
  #
  # Examples:
  # Start_gaps_size('-----seqstart--seqend') = 5
  # Start_gaps_size('-----seq') = 5
  # Start_gaps_size('seqstart----seqend') = 0
  #
  ##################
  
  dash_distribution <- grep(string[[1]], pattern = '-')
  gap_end <- c()
  if(1 %in% dash_distribution){
    for(i in c(1:(length(dash_distribution)-1))){
      if(dash_distribution[i+1] != dash_distribution[i]+1){
        gap_end <- i
        break
      }
    }
  }else{gap_end <- 0}
  if(is.null(gap_end)){gap_end <- length(dash_distribution)}
  return(gap_end)
  
}

Fix_gaps <- function(aln, pattern, subject){
  
  ##################
  # Description:
  #
  # Introduce leading and trailing gaps on alignments of class "PairwiseAlignmentsSingleSubject"
  # From Biostrings.
  #
  ##################
  
  query <- as.character(aln@pattern)
  ref <- as.character(aln@subject)
  
  # adding the leading gaps:
  
  if(start(aln@subject) != 1 | start(aln@pattern) != 1){
    if(start(aln@subject) > start(aln@pattern)){
      
      aligned_query <- as.character( aln@pattern )
      missing_query <- paste(rep(c('-'),start(subject(aln))-1),collapse = '')
      query = str_c(missing_query, aligned_query, collapse = '')
      
      aligned_ref <- as.character( aln@subject )
      missing_ref <- str_sub(as.character(subject), start = 1, end =start(subject(aln))-1)
      ref = str_c(missing_ref, aligned_ref, collapse = '')
      
    }else if(start(aln@subject) < start(aln@pattern)){
      
      aligned_ref <- as.character( aln@subject )
      missing_ref <- paste(rep(c('-'),start(pattern(aln))-1),collapse = '')
      ref = str_c(missing_ref, aligned_ref, collapse = '')
      
      aligned_query <- as.character( aln@pattern )
      missing_query <- str_sub(as.character(pattern), start = 1, end =start(pattern(aln))-1)
      query = str_c(missing_query, aligned_query, collapse = '')
      
    }
  }
  
  # adding the trailing gaps:
  
  missing_trails <- c(nchar(pattern)-end(aln@pattern), nchar(subject)-end(aln@subject))
  
  if(any(missing_trails != 0)){
    if(missing_trails[1] < missing_trails[2]){
      
      aligned_query <- query
      missing_query <- paste(rep(c('-'),missing_trails[2]),collapse = '')
      query <- str_c(aligned_query, missing_query, collapse = '')
      
      aligned_ref <- ref
      missing_ref <- str_sub(as.character(subject), start = end(subject(aln))+1, end = nchar(subject))
      ref = str_c(aligned_ref, missing_ref, collapse = '')
      
    }else if(missing_trails[1] > missing_trails[2]){
      
      aligned_ref <- ref
      missing_ref <- paste(rep(c('-'),missing_trails[1]),collapse = '')
      ref <- str_c(aligned_ref, missing_ref, collapse = '')
      
      aligned_query <- query
      missing_query <- str_sub(as.character(pattern), start = end(pattern(aln))+1, end = nchar(pattern))
      query = str_c(aligned_query, missing_query, collapse = '')
      
    }
  }
  final_query <- strsplit(query, '')
  final_ref <- strsplit(ref, '')
  
  return(list('query' = final_query,'ref' = final_ref))
  
}

Position_vector <- function(input, ref, aln){
  
  # #################
  # Description:
  #
  # Generates suitable positioning vector (1 10 20...etc) based upon a reference 
  # sequence provided.
  #
  # #################
  
  
  v_length = nchar(input)
  
  tmp = seq(10,100, 10)
  tmp = paste(tmp,collapse = '-')
  # '-' are replaced by blank spaces to avoid representation of ALL positions, 
  # just every 10 - '|' is printed every 5 positions.
  
  pos1 = str_replace_all(tmp, '-', paste("   ","|", "    ", sep = '', collapse = ''))
  
  if (v_length > 100 ){ #After position 100 number of spaces must be readjusted 
    #to 7.
    
    tmp = seq(100,1000, 10)
    tmp = paste(tmp,collapse = '-')
    pos2 = str_replace_all(tmp, '-', paste("  ","|","    ",sep = '', collapse = ''))
    
    str_sub(pos1, -3) <- '' #Position 100 is repeated by strings pos1 and pos2
    pos = str_c(pos1, pos2)
    
  }else{pos = pos1} #Ref could be lower than 100.
  
  pos = str_c(c('1   |    '), pos, sep = '') #We were missing position 1
  
  # If there is a leading sequence of '-' (gaps) we must NOT have these into 
  # account when printing position 1.
  
  lead_dashes <- paste(rep(c(' '),start(aln@pattern)-1), collapse = '')
  pos <- str_c(lead_dashes, pos, collapse = '')
  
  # v_length = lengthFixing( v_length, pos) #Length needs to be adjusted in case 
  #a number is being chopped in half. ex: length(ref)=101
  
  str_sub( pos, v_length+1, nchar( pos )) <- ''
  
  #Once positioning vector is completed '-' must be reintegrated to account for 
  #inserts with respect to ref.
  
  loc <- which( str_split(ref,pattern = '')[[1]] %in% '-')
  
  for ( i in loc ){
    pos <- unlist(strsplit(pos, split = ''))
    counter = i
    if (pos[counter] != ' '){
      if (pos[counter-1] == ' '){
        
        pos <- base::append( pos, ' ', after = counter -1)
        pos[ length( pos )] <- ''
        
      }else{
        
        while ( pos[counter] != ' '){ #-To prevent the splitting of whole numbers.
          counter <- counter + 1
        }
        #### Introduced spacers after the exit of the while loop ###
        pos <- base::append( pos, ' ', after = counter -1)
        pos[ length( pos )] <- ''
        ############################################################.
      }
    }else{
      
      pos <- base::append( pos, ' ', after = counter -1)
      pos[ length( pos )] <- ''
      
    }
  }
  
  return(pos)
}

introduce_ws_for_inserts <- function(pos_vector, aln, mut_positions){
  
  pos_vector[!str_split(aln@subject, '')[[1]] %in% '-'][mut_positions] <- mut_positions
  return(pos_vector)
  
}

lengthFixing <- function(v_length, pos){
  #############.
  # Description:
  #
  # Helper function for Position_vector()
  #
  #############.
  
  if (str_sub(pos, v_length, v_length) != ' '){
    v_length = v_length+1
    return(lengthFixing(v_length, pos))
  }
  return(v_length)
}

InDel_Extraction <- function(ALN, Abundance){
  
  Del_Sizes = as_tibble(Biostrings::deletion(ALN)) %>% select(group, width)
  In_Sizes = as_tibble(Biostrings::insertion(ALN)) %>% select(group, width)
  
  Del_Data = InDel_Sizes(Abundance, Del_Sizes)
  In_Data = InDel_Sizes(Abundance, In_Sizes)
  
  colnames(Del_Data) = c('Deletions','TotalD')
  colnames(In_Data) = c('Insertions','TotalI')
  
  return(list(Del_Data, In_Data))
  
}

InDel_Sizes <- function(Abundance, InDel_Sizes){
  
  relevant_Abundance = tibble()
  
  for( value in InDel_Sizes$group ){
    
    relevant_Abundance <- append(relevant_Abundance, Abundance[value, 2])
    
  }
  
  relevant_Abundance <- flatten_dbl(relevant_Abundance)
  relevant_Abundance <- as_tibble(relevant_Abundance) #Necesary for function
  #bind_cols() used next
  
  InDel_Data <- bind_cols(InDel_Sizes, relevant_Abundance) %>% 
    select(width, value) %>% group_by(width) %>% summarise(sum(value))
  colnames(InDel_Data) <- c('Size', 'TotalID')
  
  return(InDel_Data)
  
}

Unzip_file <- function(PATH, fileGZ_path){
  
  ##############.
  # Description:
  #
  # Decompresses .gz files into PATH and return new complete file path as a string.
  #
  ##############.
  
  file_decomposed <- str_split(fileGZ_path, '/')[[1]]
  file_name <- file_decomposed[length(file_decomposed)]
  
  if(str_sub(PATH, nchar(PATH)) != '/'){ # Make sure PATH ends with a bracket
    PATH <- str_c(PATH,'/')
  }
  
  new_file <- str_c(PATH,str_remove(file_name,pattern = '.gz'))
  
  gunzip(filename = fileGZ_path, destname = new_file, remove = FALSE) 
  
  #file and prevents compressed file deletion
  # file_path = str_sub(fileGZ_path, 1, -4) #removes last letter ## legacy code
  
  return(new_file) 
}

Unravel_Positions <- function(Insert_per_loci, Insert_data){
  
  Total_Count <- c(rep(0,nrow(Insert_per_loci)))
  
  for (i in 1:nrow(Insert_data)){
    
    start <- Insert_data['start'][i,]
    end <- Insert_data['end'][i,]
    Abundance <- Insert_data['Count'][i,]
    
    Total_Count[start:end] <- Abundance + Total_Count[start:end]
  }
  
  Insert_per_loci <- Insert_per_loci %>% add_column(Total_Count)
  
  return(Insert_per_loci)
  
}

Group_list <- function(data){
  
  if (typeof(data) == "list"){
    Groups <- c(1:nrow(data))
    for (i in 1:length(Groups)){ Groups[i] <- paste('Group ',i, sep = '') }
    
    Groups[nrow(data)] <- 'Other'
    
  }else{
    Groups <- c(1:length(data))
    for (i in 1:length(Groups)){ Groups[i] <- paste('Group ',i, sep = '') }
  }
  return(Groups)
}

Clusters_Alignments <- function(Pattern_Seq, Subject_Seq, gapOpening = 10, gapExtension = 0.5, type = "overlap" ){
  
  aln = pairwiseAlignment(Pattern_Seq, Subject_Seq, type = type, 
                          substitutionMatrix = nucleotideSubstitutionMatrix(match = 5, mismatch = -4),
                          gapOpening=gapOpening, 
                          gapExtension=gapExtension)

  #Changing deletions and insertion columns to show total number of indels 
  #instead of total number of indel bps
  number_of_dels <- map(as.list(deletion(aln)), ~length(.))
  number_of_ins <- map(as.list(insertion(aln)), ~length(.))
  
  tmp = data.frame(
    ID = names(Pattern_Seq),
    mismatch = nmismatch(aln),
    #alignment_length = nchar(aln), commented out because we wanted length of sequence not the length of the alignment
    score = score(aln),
    start = start(subject(aln)),
    end = end(subject(aln)),
    Deletions = purrr::flatten_dbl(number_of_dels),
    Deleted_bps = nindel(aln)@deletion[, 2],
    Insertions = purrr::flatten_dbl(number_of_ins),
    Inserted_bps = nindel(aln)@insertion[, 2]
  )
  
  tmp = tmp %>% separate(ID, c("ID","kk"), sep =" ") %>% select(-kk)
  
  ####################.
  #### FILL GAPS #####
  ####################.
  
  id <- str_split(names(unaligned(aln@pattern)), pattern = ' ')
  
  aligned_names <- c()
  for(i in id){aligned_names <- c(aligned_names, i[1])}
  pos_in_table <- match(aligned_names, tmp$ID)
  
  insertions <- as.data.frame(indel(aln)@insertion) %>% select(-group_name)
  deletions <- as.data.frame(indel(aln)@deletion) %>% select(-group_name)
  
  insertions_pos <- insertions %>% mutate(middle_pos = paste(start,":",end, sep = ''))
  deletions_pos <- deletions %>% mutate(middle_pos = paste(start,":",end, sep = ''))
  
  insertions_clst <- insertions_pos %>% mutate(cluster = insertions_pos$group) %>% 
    select(-group) %>% group_by(cluster) %>% summarise(insert_range = paste(middle_pos, collapse = ";"))
  deletions_clst <- deletions_pos %>% mutate(cluster = deletions_pos$group) %>% 
    select(-group) %>% group_by(cluster) %>% summarise(deletion_range = paste(middle_pos, collapse = ";"))
  
  tmp_with_rownames <- tmp %>% rownames_to_column('cluster') %>% mutate(cluster = as.numeric(cluster))
  
  tmp_with_ins <- left_join(tmp_with_rownames, insertions_clst, by = 'cluster')
  tmp_with_indels <- left_join(tmp_with_ins, deletions_clst, by = 'cluster')
  tmp <- tmp_with_indels %>% select(-cluster)
  
  # Since we have artificially included end gaps into out alignments, we have to compute the new deletions and insertions counts.
  for(i in c(1:length(aln))){
    alignment = Fix_gaps(aln[i], Pattern_Seq[i], Subject_Seq)
    end_gaps_query <- c(alignment$query[[1]][1] == '-', alignment$query[[1]][length(alignment$query[[1]])] == '-')
    # Query gaps
    if(sum(end_gaps_query) > 0){
      query_length = length(alignment$query[[1]])
      # Do we have gaps leading and trailing?
      if(sum(end_gaps_query) == 2){
        l_count= 1
        while(alignment$query[[1]][l_count] == '-'){
          l_count = l_count + 1
        }
        t_count <- query_length
        while(alignment$query[[1]][t_count] == '-'){
          t_count = t_count - 1
        }
        tmp$deletion_range[i] <- str_c('1:',l_count,';',tmp$deletion_range[i], ";",t_count,':',query_length, sep = '')
      # Only leading gaps
      }else if(end_gaps_query[1] && !end_gaps_query[2]){
        l_count = sum(alignment$query[[1]] == '-')-tmp$Deleted_bps[i]
        tmp$deletion_range[i] <- str_c('1:',l_count,';',tmp$deletion_range[i], sep = '')
      # Only trailing gaps
      }else if(!end_gaps_query[1] && end_gaps_query[2]){
        t_size = sum(alignment$query[[1]] == '-')-tmp$Deleted_bps[i]
        t_count = query_length - t_size
        tmp$deletion_range[i] <- str_c(tmp$deletion_range[i], ";",t_count,":",query_length, sep = '')
      }
      tmp$Deleted_bps[i] <- sum(alignment$query[[1]] == '-')
      tmp$Deletions[i] <- tmp$Deletions[i]+sum(end_gaps_query)
    }
    # Ref gaps
    end_gaps_ref <- c(alignment$ref[[1]][1] == '-', alignment$ref[[1]][length(alignment$ref[[1]])] == '-')
    
    if(sum(end_gaps_ref) > 0){
      ref_length = length(alignment$ref[[1]])
      #Gaps leading and trailing
      if(sum(end_gaps_ref) == 2){
        l_count = 1
        while(alignment$query[[1]][l_count] == '-'){
          l_count = l_count + 1
        }
        t_count <- ref_length
        while(alignment$ref[[1]][t_count] == '-'){
          t_count = t_count - 1
        }
        tmp$insert_range[i] <- str_c('1:',l_count,';',tmp$insert_range[i],';',t_count,':',ref_length, sep = '')
      } # Only leading gaps
    else if(end_gaps_ref[1] && !end_gaps_ref[2]){
      l_count = sum(alignment$ref[[1]] == '-')-tmp$Inserted_bps[i]
      tmp$insert_range[i] <- str_c('1:',l_count,';',tmp$insert_range[i], sep = '')
        # Only trailing gaps
    }else if(!end_gaps_ref[1] && end_gaps_ref[2]){
      t_count = sum(alignment$ref[[1]] == '-')-tmp$Inserted_bps[i]
      tmp$insert_range[i] <- str_c(tmp$insert_range[i], ";",t_count,":",ref_length, sep = '')
    }
      tmp$Inserted_bps[i] <- sum(alignment$ref[[1]] == '-')
      tmp$Insertions[i] <- tmp$Insertions[i]+sum(end_gaps_ref)
    }
  }
  
  return(list(tmp = tmp, aln = aln))
}

file.dir <- function(file_datapath){
  
  # #########.
  # Description:
  #
  # Splits and keeps the directory of a string containing a directory & a file
  #
  # ########.
  #
  # REQUIRES: ####.
  # 
  # strngr or tidyverse
  #
  # #########.
  
  str_sub(file_datapath, start = 1, end = nchar(file_datapath)-(nchar(str_split(file_datapath, pattern = '/')[[1]][length(str_split(file_datapath, pattern = '/')[[1]])])+1))
  
}

find_target_location <- function(Tabla, TablaT, Target){
  
  # ######################################
  # Description:
  #
  # Finds the target sequence using score, mismatches, indels and coverage.
  #
  # ######################################
  
  perfectMatches <- filter(TablaT, mismatch == 0, Deletions == 0, Insertions == 0, Length == nchar(Target))
  
  if(nrow(perfectMatches) == 0){
    
    Target_loc <- "Target exact match not found"
    Target_text <- " :WARNING: "
    
  }else if(nrow(perfectMatches) > 1){
    Target_loc <- NULL
    perfectMatches <- perfectMatches %>% arrange(desc(Abundance))
    
    Target_loc <- grep(Tabla$ID, pattern = head(perfectMatches)[1,1])
    Target_loc <- rownames(Tabla)[Target_loc]
    
    if(is_empty(Target_loc)){Target_loc <- 'NONE'}
    Target_text <- 'Targets found in:'
    
  }else if(nrow(perfectMatches) == 1){
    
    MatchTargetHighScore_ID <- perfectMatches$ID[1]
    Target_loc <- detect_index(Tabla$ID, ~. == MatchTargetHighScore_ID)
    Target_loc <- rownames(Tabla)[Target_loc]
    
    if(is_empty(Target_loc)){Target_loc <- 'NONE'}
    
    Target_text <- 'Target is found in:'
  }
  
  return(c(Target_text, Target_loc))
  
}

conditional_input_decision_tree <- function(path_to_files){
  
  #### DESCRIPTION: #######################################################.
  #
  # This funciton decompresses files input through detection of gz extension in combination
  # with .fastq extension. It is designed to work in the environment of mosaic finder.
  # 
  #### USAGE: #####
  #
  # ./conditional_input_decision_tree( path_to_files )
  # 
  #### REQUIRES: ####
  #
  # stringr or tidyverse
  #
  #########################################################################.
  
  file_names.split <- str_split(path_to_files$name, pattern = '\\.')[[1]] # Take one sample of the fed entries.
  path_to_files.split.extension <- file_names.split[2:length(file_names.split)] # Isolate extensions
  
  path_to_files.split <- str_split(path_to_files$datapath, pattern = '\\.')[[1]]
  path_to_files.split.path <- file.dir(path_to_files.split[1])
  
  
  ##############################################################################
  ######################## CASE .FASTQ #########################################
  ##############################################################################
  
  # If more than one file is given we now we are not being fed a zipped directory
  if(nrow(path_to_files) > 1){
    if( length(path_to_files.split.extension) == 1 && path_to_files.split.extension[1] == 'fastq'){
      
      file.rename(path_to_files$datapath, str_c(path_to_files.split.path,'/',path_to_files$name))
      new_file_paths <- dir(path_to_files.split.path, pattern = 'fastq' ,full.names = TRUE)
      
      return(new_file_paths) # There is no decompression required, return fed arguments as is.
      
      ##########################################################################
      ######################## CASE .FASTQ.GZ ##################################
      ##########################################################################
      
    }else if(length(path_to_files.split.extension) > 1 && 
             path_to_files.split.extension[length(path_to_files.split.extension)] == 'gz' && 
             path_to_files.split.extension[length(path_to_files.split.extension)-1] == 'fastq'){
      for(i in c(1:nrow(path_to_files))){
        # This decompresses files into a desired location
        gunzip(filename = path_to_files$datapath[i], destname = str_c(path_to_files.split.path,'/',str_remove(path_to_files$name[i], pattern = '.gz')), remove = FALSE)
        # We are not interested in the returned string though as we will obtain the path later.
        # file.rename(str_remove(path_to_files$datapath[i], pattern = '.gz'), path_to_files$name[i])
      }
      
      new_file_paths <- dir(path_to_files.split.path, pattern = 'fastq' ,full.names = TRUE)
      return(new_file_paths)
    }
    
    ############################################################################
    ######################## CASE .ZIP #########################################
    ############################################################################
    
  }else if(length(path_to_files.split.extension) == 1 && path_to_files.split.extension == 'zip'){
    
    if(nrow(path_to_files) > 1){return('::ERROR:: only one .zip file permitted')}
    unzip(path_to_files$datapath, exdir = path_to_files.split.path, overwrite = TRUE)
    new_file_path <- dir(path_to_files.split.path, pattern = 'fastq$',full.names = TRUE)
    
    if(is_empty(new_file_path)){
      new_file_path <- dir(path_to_files.split.path, pattern = 'fastq',full.names = TRUE)
      if(is_empty(new_file_path)){return('::ERROR:: zip file does not contain a fastq files.')}
    }
       
       if('gz' %in% str_split(new_file_path, pattern = "\\.")[[1]]){      
          for( file in new_file_path ){
            gunzip(filename = file, destname = str_c(str_remove(file, pattern = '.gz')), 
                   remove = FALSE, overwrite = TRUE) # This decompresses files into a desired location
            # We are not interested in the returned string though as we will obtain the path later.
            new_file_path <- dir(path_to_files.split.path, pattern = 'fastq$', full.names = TRUE)
            # file.rename(str_remove(path_to_files$datapath[i], pattern = '.gz'), path_to_files$name[i])
          }
    }
    
    return( new_file_path )
    
    ############################################################################
    ######################## CASE .TAR.GZ ######################################
    ############################################################################
    ##################### Should we remove this? ###############################
    ############################################################################
    
  }else if('tar' %in% path_to_files.split.extension){
    
    untar(path_to_files$datapath, exdir = path_to_files.split.path)
    new_file_path <- dir(path_to_files.split.path, pattern = 'fastq', full.names = TRUE)
    
    if(!is.na(str_split(new_file_path, pattern = "\\.")[[1]][3])){
      
      for( file in new_file_path ){
        gunzip(filename = file, destname = str_c(path_to_files.split.path,'/', str_remove(file, pattern = '.gz')), remove = FALSE) # This decompresses files into a desired location
        # We are not interested in the returned string though as we will obtain the path later.
        new_file_path <- dir(path_to_files.split.path, pattern = 'fastq$' ,full.names = TRUE)
        # file.rename(str_remove(path_to_files$datapath[i], pattern = '.gz'), path_to_files$name[i])
      }
    } # Si la extension es la correcta devuelve el string con los paths, si no que descomprima como la opción de entrada múltiple
    
    return( new_file_path )
  }
  return(str_c('Unrecognised input data: ', path_to_files))
}

read.primers <- function(path){
  
  library(tidyverse)
  library(Biostrings)
  
  F_Primer1 <- read_lines(path)[2]
  R_Primer1 <- read_lines(path)[4]
  
  F_Primer2 <- as.character(Biostrings::reverseComplement(DNAString(R_Primer1)))
  R_Primer2 <- as.character(Biostrings::reverseComplement(DNAString(F_Primer1)))
  
  store.primers <- list('F_Primer1' = F_Primer1, 'R_Primer1' = R_Primer1, 'F_Primer2' = F_Primer2, 'R_Primer2' = R_Primer2)
  return(store.primers)
}

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
