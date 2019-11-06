scQC <- function(X){
  X[X < 0] <- 0
  lSize <- apply(X,2,sum)
  X <- X[,lSize > 1000]
  lSize <- apply(X,2,sum)
  X <- X[,!lSize %in% boxplot.stats(lSize)$out]
  mtRate <- X[grepl('MT-',toupper(rownames(X))),]
  mtRate <- apply(mtRate,2,sum)/apply(X,2,sum)
  X <- X[,mtRate < 0.1]
  X <- as(X, 'dgCMatrix')
  return(X)
}
