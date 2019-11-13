#' @export cpmNormalization
#' @title cpmNormalization
#' @description Perfrom counts per million (CPM) normalization
#' @param X ...
#' @return ...
cpmNormalization <- function(X){
  X <- as.matrix(X)
  X <- t(t(X)/colSums(X))*1e6
  X <- as(X, 'dgCMatrix')
  return(X)
}
