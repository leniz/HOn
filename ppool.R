p.pool <-function(pValues,u=1) { 
  # this functions pools the pValues using the partial conjunction p-value based on Simes
  # pValues is a vector
  m = length(pValues)
  pV = sort(pValues)
  pV = pV[u:m]
  m = m - u + 1
  pV = m*pV/c(1:m)
  
  return(min(pV))
}
