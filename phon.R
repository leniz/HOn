phon <- function(pValues,
                     alpha = 0.05,
                     kappa = 0.001){
  ## pool the p-values, partial conjunction
  depth <- list.depth(pValues)
  pconj <- vector("list",depth)
  pconj.aux <- pValues
  ml <- vector("list", depth)
  for (j in 1:(depth-1)){
   ml[[j]] <- unlist(lapply(pconj.aux,length))
   pconj.aux <- unlist(pconj.aux, recursive = F)
  }
  
  ml[[depth]] <- unlist(lapply(pconj.aux,length)) # size of the families
  stru <- lapply(ml, function(x) rep(1:length(x),x))
  # obtain u_l for partial conjunction
  ul <- as.integer(ml[[depth]]*kappa)
  ul[ul<ml[[depth]]*kappa] <- ul[ul<ml[[depth]]*kappa] + 1
  # pool pvalues of families
  for (j in 1:length(pconj.aux)){pconj[[depth]] <- c(pconj[[depth]],p.pool(pconj.aux[[j]], u=ul[j]) )}
  # pool pooled pvalues
  for (jj in 1:(depth-1)){
    idx <- depth-jj
    for (j in 1:length(ml[[idx]])){pconj[[idx]] <- c(pconj[[idx]], p.pool(pconj[[idx+1]][stru[[idx]]==j], u=1) )}
  }
  
  ## screen 
  screen <- lapply(pconj, function(x) x<=alpha*kappa)
  idx <- screen[[1]]
  for (j in 2:depth){
    idx<-rep(idx,ml[[j-1]])
    idx<-idx*screen[[j]]
  }
  
return(list(ml, pconj.aux,ul,pconj,screen,idx))
}
