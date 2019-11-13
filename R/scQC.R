#' @export scQC
#' @title scQC
#' @importFrom grDevices boxplot.stats
#' @description ...
#' @param X ...
#' @param minLibSize ...
#' @param removeOutlierCells ...
#' @param minPCT ...
#' @param maxMTratio ...
#' @return ...
#' @references ...
scQC <- function(X, minLibSize = 1000, removeOutlierCells = TRUE, minPCT = 0.05, maxMTratio = 0.1){
  # Removing values lower than 0
  X[X < 0] <- 0
  # Computing library size
  lSize <- apply(X,2,sum)
  # Filtering out by minimum library size
  X <- X[,lSize > minLibSize]
  # Removing outlier cells
  if(removeOutlierCells){
    lSize <- apply(X,2,sum)
    X <- X[,!lSize %in% boxplot.stats(lSize)$out]
  }
  # Computing mitochondrial ratio
  mtRate <- X[grepl('MT-',toupper(rownames(X))),]
  mtRate <- apply(mtRate,2,sum)/apply(X,2,sum)
  # Filtering out by mitochondrial ratio
  X <- X[,mtRate < maxMTratio]
  # Filtering out by minPCT
  X <- X[apply(X!=0,1,mean) > minPCT,]
  # Setting the output structure
  X <- as(X, 'dgCMatrix')
  # Return
  return(X)
}
