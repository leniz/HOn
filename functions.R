# this functions need the following functions: 
# source('~/R/codeANDdata/hierBH.R')
# source('~/R/codeANDdata/pVpool.R')
# source('~/R/codeANDdata/hierasymptkappa.R')


# Simes pooled p-values--------------------------------------------
# pools pValues using Simes 

p.pool <-function(pValues,u=1) { 
  # later: add methods
  # pValues is a vector
  m = length(pValues)
  pV = sort(pValues)
  pV = pV[u:m]
  m = m - u + 1
  pV = m*pV/c(1:m)
  
  return(min(pV))
}



# Hierarchically ordered procedure for pvalues without family structure ------------------------------

phiHO_1 <- function(pValues, # 1 x m vector of pValues
                    alpha = 0.05,
                    kappa = 0.001) {
  m <- length(pValues)
  u <- as.integer(kappa*m) + 1
  fpvalue <- p.pool(pValues, u) # family p-value
  if (fpvalue <= alpha*kappa) { # pretest of the one family p-value
    hresult <- aorc(pValues, alpha, startIDX_SUD = u, betaAdjustment = 0, silent = TRUE) #local result in the family by AORC test
    result <- list(hresult)
  }
  else { # no rejections made
    result <- list(pValues != pValues)
  }
  return(result)
} 
# phiHO (end)



# hierarchically ordered procedure 2 stages ----------------------------- 
# as in  PLOS ONE's ::KS,KT,TD 2016:: for pValues with a family structure
# two stages
phiHO_2 <- function(pValues, 
                    alpha = 0.05, 
                    kappa = 0.001){
  result <- list()
  if (is.list(pValues) == FALSE) { # no family structure, test as one single family
    result <- phiHO_1(pValues, alpha, kappa)
  } else {
    k <- length(pValues)  # number of families
    if (1/k <= kappa) kappa <- 1/k # For larger families obtain a smaller kappa than the preset
    for (i in 1:k) {
      if (length(pValues[[i]]) == 0) {
        lresult <- list(FALSE)
        result <- c(result, lresult)
        next
      }
      fpValues <- pValues[[i]] # family p-values
      fm <- length(fpValues)
      fu <- as.integer(kappa*fm) # determine u_i
      if (fu < kappa*fm) fu <- fu + 1
      fpvalue <- p.pool(fpValues, fu) # family p-value
      if (fpvalue <= alpha*kappa) { # pretest of the family p-value
        lresult <- aorc(fpValues, alpha, startIDX_SUD = fu, betaAdjustment = 0, silent = TRUE) # local test of the family by AORC SUD procedure
        lresult <- list(lresult$rejected)
      } else { # no rejection in this family is made
        lresult <- list(fpValues!=fpValues)
      } 
      result <- c(result, lresult)
    }
    return(result)
  }
}
# phi^HO2 (end)

# hierarchically ordered procedure three stages --------------
# 
phiHO_3 <- function(pValues,
                    alpha = 0.05,
                    kappa = 0.001){
  
  aux_kappa <- kappa
  result <- list()
  pool <- vector("list",length(pValues))
  
  if (!is.list(pValues)) {
    result <- phiHO_2(pValues)
  } else {
    
    # pool p-values from families of indiv. hyp (Partial conjunction hypothesis)
    for (i in 1:length(pValues)){
      x <- numeric()
      for (j in 1:length(pValues[[i]])){
        fm <- length(pValues[[i]][[j]])
        fu <- as.integer(kappa*fm) # determine u_i
        if (fu < kappa*fm) fu <- fu + 1
        fu <- 1 # SOLO PARA PRUEBAS. QUITAR PRONTO
        x <- c(x,p.pool(pValues[[i]][[j]], fu))
      }
      pool[[i]] <- x
    }
    
    
    # select superfamilies to test pooling again the above p-values
    k <- length(pValues)
    if (1/k <= kappa) kappa <- 1/k # For larger families obtain a smaller kappa than the preset
    for (i in 1:k) {
      if (length(pValues[[i]]) == 0) {
        lresult <- list(FALSE)
        result <- c(result, lresult)
        next
      }
      
      # we use for now simes with u=1, because we consider there will be a small number 
      # of families of indiv. p-values in each of the upper level families
      # fu <- as.integer(kappa*fm) # determine u_i
      # if (fu < kappa*fm) fu <- fu + 1
      fu <- 1
      fpvalue <- p.pool(pool[[i]],u=fu)
      if (fpvalue <= alpha*kappa) { # pretest of the family p-value (fpvalue <= alpha*kappa)
        lresult <- phiHO_2(pValues[[i]],alpha, aux_kappa)
      } else { # no rejection in this family is made
        lresult <- list(fpvalue!=fpvalue)
      }
      result <- c(result, list(lresult))
      
    }
    
  }
  return(result)
}




list.depth <- function(this, thisdepth = 0) {
  # http://stackoverflow.com/a/13433689/1270695
  if(!is.list(this)) {
    return(thisdepth)
  } else {
    return(max(unlist(lapply(this, list.depth, thisdepth = thisdepth+1))))    
  }
}
