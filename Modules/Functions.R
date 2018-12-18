

Position_vector <- function(input, ref, primers = TRUE){
  #Generates suitable positioning vector (1 10 20...etc) based upon a reference 
  #sequence provided.
  
  v_length = nchar(input)
  
  tmp = seq(10,100, 10)
  tmp = paste(tmp,collapse = '-')
  # '-' are replaced by blanck spaces to avoid representation of ALL position, 
  # just every 10
  pos1 = str_replace_all(tmp, '-', paste(rep(x = ' ',8), collapse = ''))

  if (v_length > 100 ){ #After position 100 number of spaces must be readjusted 
    #to 7.
    
    tmp = seq(100,1000, 10)
    tmp = paste(tmp,collapse = '-')
    pos2 = str_replace_all(tmp, '-', paste(rep(x = ' ',7), collapse = ''))
    
    str_sub(pos1, -3) <- '' #Position 100 is repeated by strings pos1 and pos2
    pos = str_c(pos1, pos2)
    
  }else{pos = pos1} #Ref could be lower than 100.
  
  pos = str_c(c('1        '), pos) #We were missing position 1
  
  v_length = lengthFixing( v_length, pos) #Length needs to be adjusted in case 
  #a number is being chopped in half. ex: length(ref)=101
  
  str_sub( pos, v_length+1, nchar( pos )) <- ''
  
  #Once positioning vector is completed '-' must be reintegrated to account for 
  #inserts with respect to ref.
  
  loc <- which( ref %in% '-')

  for ( i in loc ){
    pos <- unlist(strsplit(pos, split = ''))
    counter = i
    while ( pos[counter] != ' '){ #-To prevent the splitting of whole numbers.
      counter <- counter + 1
    }
    pos <- base::append( pos, ' ', after = counter-1 )
    pos[ length( pos )] <- ''
    #pos[length(pos)] <- ''
}
  
  return(pos)
}

lengthFixing <- function(v_length, pos){
  if (str_sub(pos, v_length, v_length) != ' '){
    v_length = v_length+1
    return(lengthFixing(v_length, pos))
  }
  return(v_length)
}

InDel_Extraction <- function(ALN, Abundance){
  
  Del_Sizes = as.data.frame(nindel(ALN)@deletion[ ,2])
  In_Sizes = as.data.frame(nindel(ALN)@insertion[ ,2])
  
  Del_Data = InDel_Sizes(Abundance, Del_Sizes)
  In_Data = InDel_Sizes(Abundance, In_Sizes)
  
  colnames(Del_Data) = c('TotalD','Deletions')
  colnames(In_Data) = c('TotalI','Insertions')
  
  return(list(Del_Data, In_Data))
  
}

InDel_Sizes <- function(Abundance, InDel_Sizes){
  
  InDel_Sizes = cbind.data.frame(rep.int(1, nrow(InDel_Sizes)), InDel_Sizes)
  colnames(InDel_Sizes) = c('Amount', 'InDels')
  
  InDel_Data = bind_cols(InDel_Sizes, Abundance) %>% 
    mutate(Total = Abundance*Amount) %>% 
    select(Total, InDels) %>% group_by(InDels) %>% 
    mutate_at(vars(Total), funs(sum(Total))) %>% 
    group_by(Total) %>% distinct(InDels) %>% filter(InDels > 0)
  
  return(InDel_Data)
  
}

Unzip_file <- function(PATH, fileGZ_path){
  #Decompresses .gz files and return new file path as string.
  
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

Clusters_Alignments <- function(Clust_Seq, Align_Seq){
  
  aln = pairwiseAlignment(Clust_Seq, Align_Seq, type = "global", gapOpening=10, 
                          gapExtension=5)
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