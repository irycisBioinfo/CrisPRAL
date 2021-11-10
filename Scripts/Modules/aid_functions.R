#-Contains several R functions with various porpouses
library(tidyverse)

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

CM_generator <- function(list1, list2){
#Make a contact matrix form two list of strings.

  emptyM <- matrix(nrow = nrow(list1), ncol = nrow(list2))
  emptyM[,] <- 0
  rownames(emptyM) <- list1[[1]]
  colnames(emptyM) <- list2[[1]]
  
  

    for(entry in intersect(list1[[1]], list2[[1]])){
      
      # print(c(match(gene, rownames(emptyM)),match(gene, colnames(emptyM))))
      emptyM[match(entry, rownames(emptyM)),match(entry, colnames(emptyM))] <- 1
      
      
    }
  return(emptyM)
}

dir.create_force <- function(dir.name, iter = 0){
  while(dir.exists(dir.name)){
    iter <- iter + 1
    if(iter == 1){
      dir.name <- str_c(dir.name, "_", iter)}
    else{dir.name <- str_replace(dir.name, paste('_', as.character(iter-1), sep = ''), 
                                 paste('_',as.character(iter), sep = ''))}
  }
  
  dir.create(dir.name)
  return(dir.name)}

intersect_4 <- function(data_list = NULL, data1 = NULL, data2 = NULL, data3 = NULL, data4 = NULL){
  result <- "No data was entered, pls provide a list with 4 data sets"
  if(!is.null(data1)){
    result <- dplyr::intersect(dplyr::intersect(dplyr::intersect(data1, data2),data3), data4)
  }
  
  if(!is.null(data_list)){
    result <- dplyr::intersect(dplyr::intersect(dplyr::intersect(data_list[[1]], data_list[[2]]),data_list[[3]]), data_list[[4]])
  }
  
  return(result)
}

intersect_3 <- function(data_list = NULL, data1 = NULL, data2 = NULL, data3 = NULL){
  result <- "No data was entered, pls provide a list with 4 data sets"
  if(!is.null(data1)){
    result <- dplyr::intersect(dplyr::intersect(dplyr::intersect(data1, data2),data3), data4)
  }
  
  if(!is.null(data_list)){
    result <- dplyr::intersect(dplyr::intersect(data_list[[1]], data_list[[2]]),data_list[[3]])
  }
  
  return(result)
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

str_replace_last <- function(string, pattern = '', replacement = ''){
  
  ########################
  # Description:
  #
  # This funciton allows the use str_replace() from tidyverse to replace only the
  # LAST match found, opposite to its main function that only replaces the first
  # match.
  #
  # Requirements:
  #  
  # library(tidyverse)
  #
  ########################
  
  if(is_empty(pattern)){return(print('Pattern may not be empty'))}
  
  return(reverse(str_replace(reverse(string), reverse(pattern), reverse(replacement))))
  
}
