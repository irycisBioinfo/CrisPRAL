

Position_vector <- function(input, ref, primers = TRUE){
  #Generates suitable positioning vector (1 10 20...etc) based upon a reference sequence provided.
  
  v_length = nchar(input)
  
  tmp = seq(10,100, 10)
  tmp = paste(tmp,collapse = '-')
  # '-' are replaced by blanck spaces to avoid representation of ALL position, just every 10
  pos1 = str_replace_all(tmp, '-', paste(rep(x = ' ',8), collapse = ''))

  if (v_length > 100 ){ #After position 100 number of spaces must be readjusted to 7.
    
    tmp = seq(100,1000, 10)
    tmp = paste(tmp,collapse = '-')
    pos2 = str_replace_all(tmp, '-', paste(rep(x = ' ',7), collapse = ''))
    
    str_sub(pos1, -3) <- '' #Position 100 is repeated by strings pos1 and pos2
    pos = str_c(pos1, pos2)
    
  }else{pos = pos1} #Ref could be lower than 100.
  
  pos = str_c(c('1        '), pos) #We were missing position 1
  
  v_length = lengthFixing( v_length, pos) #Length needs to be adjusted in case a number is being chopped in half. ex: length(ref)=101
  
  str_sub( pos, v_length, nchar( pos )) <- ''
  
  #Once positioning vector is completed '-' must be reintegrated to account for inserts with respect to ref.
  
  loc <- which( ref %in% '-')

  for (i in loc){
    pos <- unlist(strsplit(pos, split = ''))
    pos <- base::append( pos, ' ', after = i-1 )
    pos[length(pos)] <- ''
}
  
  return(pos)
}

lengthFixing <- function(v_length, pos){
  if (str_sub(pos, v_length, v_length) != ' '){
    v_length = v_length+1
    return(lengthFixing(v_length, pos))
  }
  return(v_length+1)
}