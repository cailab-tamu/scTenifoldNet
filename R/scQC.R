#' @export scQC
#' @title Single-cell quality control
#' @importFrom grDevices boxplot.stats
#' @importFrom Matrix colSums rowMeans
#' @description This function performs quality control filters over the provided input matrix, the function checks for minimum cell library size, mitochondrial ratio, outlier cells, and the fraction of cells where a gene is expressed.
#' @param X Raw counts matrix with cells as columns and genes (symbols) as rows.
#' @param minLibSize An integer value. Defines the minimum library size required for a cell to be included in the analysis.
#' @param removeOutlierCells A boolean value (TRUE/FALSE), if TRUE, the identified cells with library size greater than 1.58 IQR/sqrt(n) computed from the sample, are removed. For further details see: \code{?boxplot.stats}
#' @param minPCT A decimal value between 0 and 1. Defines the minimum fraction of cells where the gene needs to be expressed to be included in the analysis.
#' @param maxMTratio A decimal value between 0 and 1. Defines the maximum ratio of mitochondrial reads (mithocondrial reads / library size) present in a cell to be included in the analysis. It's computed using the symbol genes starting with 'MT-' non-case sensitive. 
#' @return A dgCMatrix object with the cells and the genes that pass the quality control filters.
#' @references Ilicic, Tomislav, et al. "Classification of low quality cells from single-cell RNA-seq data." Genome biology 17.1 (2016): 29.
#' @examples 
#' library(scTenifoldNet)
#' 
#' # Simulating of a dataset following a negative binomial distribution with high sparcity (~67%)
#' nCells = 2000
#' nGenes = 100
#' set.seed(1)
#' X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
#' X <- round(X)
#' X <- matrix(X, ncol = nCells)
#' rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
#' 
#' # Performing Single cell quality control
#' qcOutput <- scQC(
#'   X = X,
#'   minLibSize = 30,
#'   removeOutlierCells = TRUE,
#'   minPCT = 0.05,
#'   maxMTratio = 0.1
#' )
#' 
#' # Visualizing the Differences
#' par(
#'   mfrow = c(2, 2),
#'   mar = c(3, 3, 1, 1),
#'   mgp = c(1.5, 0.5, 0)
#' )
#' # Library Size
#' plot(
#'   Matrix::colSums(X),
#'   ylim = c(20, 70),
#'   ylab = 'Library Size',
#'  xlab = 'Cell',
#'   main = 'Library Size - Before QC'
#' )
#' abline(h = c(30, 58),
#'        lty = 2,
#'        col = 'red')
#' plot(
#'   Matrix::colSums(qcOutput),
#'   ylim = c(20, 70),
#'   ylab = 'Library Size',
#'   xlab = 'Cell',
#'   main = 'Library Size - After QC'
#' )
#' abline(h = c(30, 58),
#'        lty = 2,
#'        col = 'red')
#' # Mitochondrial ratio
#' mtGenes <- grepl('^mt-', rownames(X), ignore.case = TRUE)
#' plot(
#'   Matrix::colSums(X[mtGenes,]) / Matrix::colSums(X),
#'   ylim = c(0, 0.3),
#'   ylab = 'Mitochondrial Ratio',
#'   xlab = 'Cell',
#'   main = 'Mitochondrial Ratio - Before QC'
#' )
#' abline(h = c(0.1), lty = 2, col = 'red')
#' plot(
#'   Matrix::colSums(qcOutput[mtGenes,]) / Matrix::colSums(qcOutput),
#'   ylim = c(0, 0.3),
#'   ylab = 'Mitochondrial Ratio',
#'   xlab = 'Cell',
#'   main = 'Mitochondrial Ratio - Before QC'
#' )
#' abline(h = c(0.1), lty = 2, col = 'red')

scQC <- function(X, minLibSize = 1000, removeOutlierCells = TRUE, minPCT = 0.05, maxMTratio = 0.1){
  # Removing values lower than 0
  X[X < 0] <- 0
  # Computing library size
  lSize <- Matrix::colSums(X)
  # Filtering out by minimum library size
  X <- X[,lSize > minLibSize]
  # Removing outlier cells
  if(removeOutlierCells){
    lSize <- Matrix::colSums(X)
    X <- X[,!lSize %in% boxplot.stats(lSize)$out]
  }
  # Computing mitochondrial ratio
  mtGenes <- grepl('MT-',toupper(rownames(X)))
  if(sum(mtGenes) > 0){
    mtRate <- X[mtGenes,]
    mtRate <- Matrix::colSums(mtRate)/Matrix::colSums(X)
    # Filtering out by mitochondrial ratio
    X <- X[,mtRate < maxMTratio]
  } else {
    warning('Mitochondrial genes were not found. Be aware that apoptotic cells may be present in your sample.')
  }
  # Filtering out by minPCT
  X <- X[Matrix::rowMeans(X!=0) > minPCT,]
  # Extracting gene and cell ids
  gNames <- rownames(X)
  cNames <- colnames(X)
  # Setting the output structure
  X <- as(X, 'dgCMatrix')
  # Assigning ids
  rownames(X) <- gNames
  colnames(X) <- cNames
  # Return
  return(X)
}
