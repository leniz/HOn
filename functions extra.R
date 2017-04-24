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


pr3.aux <- pr3
unl <- function(x)
  {
  x.aux <- x
  xstr <- vector("list", list.depth(x)-1)
  for (j in 1:(list.depth(x)-1)){
    xstr[[j]] <- lapply(x.aux,length)
    x.aux <- unlist(x.aux, recursive = F)
    
  }
  return(list(vecs=x.aux, str = xstr))
} 


unl(pr3)
