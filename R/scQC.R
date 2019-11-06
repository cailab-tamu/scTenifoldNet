scQC <- function(X){
  X[X < 0] <- 0
  lSize <- colSums(X)
  X <- X[,lSize > 1000]
  lSize <- colSums(X)
  X <- X[,!lSize %in% boxplot.stats(lSize)$out]
  mtRate <- X[grepl('MT-',toupper(rownames(X))),]
  mtRate <- colSums(mtRate)/colSums(X)
  X <- X[,mtRate < 0.1]
  X <- as(X, 'dgCMatrix')
  return(X)
}
