CPM <- function(X){
  X <- as.matrix(X)
  t(t(X)/colSums(X))*1e6
}
