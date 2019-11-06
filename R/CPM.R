CPM <- function(X){
  X <- as.matrix(X)
  X <- t(t(X)/colSums(X))*1e6
  X <- as(X, 'dgCMatrix')
  return(X)
}
