# define fitting function for partial linear model with 1 monotone-related variable

part_fit <- function(x, y, w, ...){
  
  # assume that monotone variable is first column in x, unless specified otherwise
  dotarg = list(...)
  if("mon_index" %in% names(dotarg)){
    ind <- mon_index
  } 
  else{
    ind <- 1
  }
  
  
  
}