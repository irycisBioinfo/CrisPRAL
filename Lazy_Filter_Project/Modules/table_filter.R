
table_filter <- function(input, output, session){
  
  for (i in 1:length(unlist(panel))){
  
    dt_temp = filter(dt, str_detect(dt[,'RegionID'],as.character(panel[i,])))
  
    if (i == 1){ dt_filtered = dt_temp
    
    }else{
      
      dt_filtered = bind_rows(dt_filtered, dt_temp)
      
    }
  }
  
  return(dt_filtered)
  
}
