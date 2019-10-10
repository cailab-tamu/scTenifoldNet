avgNetworks <- function(netList){
  X <- netList[[1]]
  for(i in seq_along(netList)[-1]){
    X <- X + netList[[i]]
  }
  X <- X/length(netList)
  return(X)
}
