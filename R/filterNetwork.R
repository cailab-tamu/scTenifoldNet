filterNetwork <- function(network, q = 0.95){
  network <- as.matrix(network)
  absNet <- abs(network)
  network[absNet < quantile(absNet, q)] <- 0
  return(network)
}
