phon <- function(pValues,
                     alpha = 0.05,
                     kappa = 0.001){
  # pool the p-values
  depth <- list.depth(pValues)
  pconj <- vector("list",length(pValues))
  pconj.aux <- pValues
  ml <- vector("list", depth)
  for (j in 1:(depth-1)){
   ml[[j]] <- unlist(lapply(pconj.aux,length))
   pconj.aux <- unlist(pconj.aux, recursive = F)
  }
  
  ml[[depth]] <- unlist(lapply(pconj.aux,length)) # size of the families

return(list(ml, pconj.aux))
}

