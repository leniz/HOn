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
  # p.pool
  for (j in 1:length(pconj.aux)){pconj[[depth]] <- p.pool(pconj.aux[[j]], u=ul[j])}
  for (jj in 1:(depth-1)){
    for (j in 1:length(ml[[2]])){pconj[[2]] <- c(pconj[[2]], p.pool(pconj[[3]][stru[[2]]==j], u=1) )}
    for (j in 1:length(ml[[1]])){pconj[[1]] <- c(pconj[[1]], p.pool(pconj[[2]][stru[[1]]==j], u=1) )}
  }
  
  # 
return(list(ml, pconj.aux,ul,pconj))
}
