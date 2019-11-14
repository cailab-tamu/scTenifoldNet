#' @export scQC
#' @title Perform single-cell quality control
#' @importFrom grDevices boxplot.stats
#' @description This function performs quality control filters over the provided input matrices, the function checks for minimum cell library size, mitochondrial ratio, outlier cells, and the fraction of cells where the gene is expressed.
#' @param X Raw counts matrix with cells as columns and genes (symbols) as rows.
#' @param minLibSize An integer value. Defines the minimum library size required for a cell to be included in the analysis.
#' @param removeOutlierCells A boolean value (TRUE/FALSE), if TRUE, the identified cells with library size greater than 1.58 IQR/sqrt(n) computed from the sample, are removed. For further details see: ?boxplot.stats
#' @param minPCT A decimal value between 0 and 1. Defines the minimum fraction of cells where the gene needs to be expressed to be included in the analysis.
#' @param maxMTratio A decimal value between 0 and 1. Defines the maximum ratio of mitochondrial reads (mithocondrial reads / library size) present in a cell to be included in the analysis. It's computed using the symbol genes starting with 'MT-'
#' @return A dgCMatrix object with the cells and the genes that pass the quality control filters.

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
