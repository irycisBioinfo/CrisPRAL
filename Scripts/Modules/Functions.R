###################
# DESCRIPTION:
#
# This script stores different multiporpouse functions required
# by different scripts. A function is declared here if the extension of it
# is larger than 5 code lines or if it serves multiple scripts.
#
# Functions:
#
# Start_gaps_size() - used in serverAlignment.R
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
# Clusters_Alignments() - used in serverAlignment.R
#
# file.dir() - used in serverMosaic.R and others
#
###################

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

Position_vector <- function(input, ref){
  
  ##################
  # Description:
  #
  # Generates suitable positioning vector (1 10 20...etc) based upon a reference 
  # sequence provided.
  #
  ##################
  
  v_length = nchar(input)
  
  tmp = seq(10,100, 10)
  tmp = paste(tmp,collapse = '-')
  # '-' are replaced by blank spaces to avoid representation of ALL positions, 
  # just every 10
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
  
  # v_length = lengthFixing( v_length, pos) #Length needs to be adjusted in case 
  #a number is being chopped in half. ex: length(ref)=101
  
  str_sub( pos, v_length+1, nchar( pos )) <- ''
  
  #Once positioning vector is completed '-' must be reintegrated to account for 
  #inserts with respect to ref.
  
  loc <- which( ref %in% '-')
  
  for ( i in loc ){
    pos <- unlist(strsplit(pos, split = ''))
    counter = i
    if (pos[counter] != ' '){
      if (pos[counter-1] == ' '){
        
        pos <- base::append( pos, ' ', after = counter -1) #Why counter -1? lets try to remove it
        pos[ length( pos )] <- ''
        
      }else{
        
        while ( pos[counter] != ' '){ #-To prevent the splitting of whole numbers.
          counter <- counter + 1
          
        }
      }
    }else{
      
      pos <- base::append( pos, ' ', after = counter -1) #Why counter -1? lets try to remove it
      pos[ length( pos )] <- ''
      
    }
  }
  return(pos)
}

lengthFixing <- function(v_length, pos){
  #############
  # Description:
  #
  # Helper function for Position_vector()
  #
  #############
 
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
  
  ##############
  # Description:
  #
  # Decompresses .gz files and return new file path as string.
  #
  ##############
  
  system(paste(PATH, '/bin/gunzip -k ', fileGZ_path, sep = '')) #decompresses 
  #file and prevents compressed file deletion
  file_path = str_sub(fileGZ_path, 1, -4) #removes last letter
  
 return(file_path) 
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

Clusters_Alignments <- function(Clust_Seq, Align_Seq, gapOpening = 10, gapExtension = 0.5 ){
  
  aln = pairwiseAlignment(Clust_Seq, Align_Seq, type = "global", 
                          substitutionMatrix = nucleotideSubstitutionMatrix(match = 5, mismatch = -4),
                          gapOpening=gapOpening, 
                          gapExtension=gapExtension)
    tmp = data.frame(
      ID = names(Clust_Seq),
      mismatch = nmismatch(aln),
      length = nchar(aln),
      score = score(aln),
      width = width(pattern(aln)),
      start = start(pattern(aln)),
      end = end(pattern(aln)),
      deletions = nindel(aln)@deletion[, 2],
      insertions = nindel(aln)@insertion[, 2]
  )
    tmp = tmp %>% separate(ID, c("ID","kk"), sep =" ") %>% select(-kk)
    
    return(list(tmp = tmp, aln = aln))
}

file.dir <- function(file_datapath){
  ##########
  # Description:
  #
  # Splits and keeps the directory of a string containing a directory & a file
  ##########
  
 str_sub(file_datapath, start = 1, end = nchar(file_datapath)-(nchar(str_split(file_datapath, pattern = '/')[[1]][length(str_split(file_datapath, pattern = '/')[[1]])])+1))
 
}

find_target_location <- function(Tabla, TablaT){
 
 perfectMatches <- filter(TablaT, mismatch == 0, Deletions == 0, Insertions == 0)
 if(nrow(perfectMatches) == 0){

 Target_loc <- "Target exact match not found"
 Target_text <- " :WARNING: "
 
}else if(nrow(perfectMatches) > 1){
 Target_loc <- NULL
 perfectMatches <- perfectMatches %>% arrange(desc(Abundance))
 
 Target_loc <- grep(Tabla$ID, pattern = head(perfectMatches)[1,1])
 
 if(is_empty(Target_loc)){Target_loc <- 'NONE'}
 Target_text <- 'Targets found in Clusters:'
 
}else if(nrow(perfectMatches) == 1){
 
 MatchTargetHighScore_ID <- perfectMatches$ID[1]
 Target_loc <- detect_index(Tabla$ID, ~. == MatchTargetHighScore_ID)
 if(Target_loc == 0){Target_loc <- 'NONE'}
 
 Target_text <- 'Target is found in Cluster:'
}
 
 return(c(Target_loc, Target_text))
 
 }