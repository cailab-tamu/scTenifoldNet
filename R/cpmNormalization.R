#' @export cpmNormalization
#' @title Perform counts per million (CPM) normalization
#' @description Each gene count for each cell are divided by the total counts for that cell and multiplied by 1e6. No log-transformation is applied.
#' @param X Raw counts matrix with cells as columns and genes (symbols) as rows
#' @return A dgCMatrix object with the count per million (CPM) normalized values.

cpmNormalization <- function(X){
  X <- as.matrix(X)
  X <- t(t(X)/colSums(X))*1e6
  X <- as(X, 'dgCMatrix')
  return(X)
}
