#' @export scTenifoldNet
#' @title scTenifoldNet
#' @importFrom methods as
#' @importFrom rTensor as.tensor cp
#' @description Gonstruct and compare single-cell gene regulatory networks (scGRNs) using single-cell RNA-seq (scRNA-seq) data sets collected from different conditions based on principal component regression, tensor decomposition, and manifold alignment.
#' @param X Raw counts matrix with cells as columns and genes (symbols) as rows.
#' @param Y Raw counts matrix with cells as columns and genes (symbols) as rows.
#' @param qc_minLibSize An integer value. Defines the minimum library size required for a cell to be included in the analysis.
#' @param qc_removeOutlierCells A boolean value (TRUE/FALSE), if TRUE, the identified cells with library size greater than 1.58 IQR/sqrt(n) computed from the sample, are removed. For further details see: \code{?boxplot.stats}
#' @param qc_minPCT A decimal value between 0 and 1. Defines the minimum fraction of cells where the gene needs to be expressed to be included in the analysis.
#' @param qc_maxMTratio A decimal value between 0 and 1. Defines the maximum ratio of mitochondrial reads (mithocondrial reads / library size) present in a cell to be included in the analysis. It's computed using the symbol genes starting with 'MT-' non-case sensitive.
#' @param nc_nNet An integer value. The number of networks based on principal components regression to generate.
#' @param nc_nCells An integer value. The number of cells to subsample each time to generate a network.
#' @param nc_nComp An integer value. The number of principal components in PCA to generate the networks. Should be greater than 2 and lower than the total number of genes.
#' @param nc_symmetric A boolean value (TRUE/FALSE), if TRUE, the weights matrix returned will be symmetric.
#' @param nc_scaleScores A boolean value (TRUE/FALSE), if TRUE, the weights will be normalized such that the maximum absolute value is 1.
#' @param nc_q A decimal value between 0 and 1. Defines the cut-off threshold of top q\% relationships to be returned.
#' @param td_K An integer value. Defines the number of rank-one tensors used to approximate the data using CANDECOMP/PARAFAC (CP) Tensor Decomposition. 
#' @param td_maxIter An integer value. Defines the maximum number of iterations if error stay above \code{td_maxError}.
#' @param td_maxError A decimal value between 0 and 1. Defines the relative Frobenius norm error tolerance.
#' @param ma_nDim An integer value. Defines the number of dimensions of the low-dimensional feature space to be returned from the non-linear manifold alignment.
#' @param dc_minDist A decimal value. Defines the cut-off threshold of distance to limit the testing to genes that show, at least \code{dc_minDist} deviation.
#' @return A list with 3 slots as follows: 
#' \itemize{
#' \item{tensorNetworks:} The generated weight-averaged denoised gene regulatory networks using CANDECOMP/PARAFAC (CP) Tensor Decomposition.
#' \item{manifoldAlignment:} The generated low-dimensional features result of the non-linear manifold alignment.
#' \item{diffCoexpression:} The results of the differential coexpression analysis.
#' }
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
#' # Generating a perturbed network modifying the expression of genes 10, 2 and 3
#' Y <- X
#' Y[10,] <- Y[50,]
#' Y[2,] <- Y[11,]
#' Y[3,] <- Y[5,]
#' 
#' \dontrun{
#' # scTenifoldNet
#' Output <- scTenifoldNet(X = X, Y = Y,
#'                        nc_nNet = 10, nc_nCells = 500,
#'                        td_K = 3, qc_minLibSize = 30,
#'                        dc_minDist = 0)
#' 
#' # Structure of the output
#' str(Output)
#' 
#' # Differential coexpression
#' head(Output$diffCoexpression,n = 10)
#' 
#' # Plotting
#' # Genes with FDR < 0.1 are labeled as red
#' geneColor <- ifelse(Output$diffCoexpression$p.adj < 0.1, 'red', 'black')
#' qqnorm(Output$diffCoexpression$Z, pch = 16, main = 'Standardized distance', col = geneColor)
#' qqline(Output$diffCoexpression$Z)}

scTenifoldNet <- function(X, Y, qc_minLibSize = 1000, qc_removeOutlierCells = TRUE,
                          qc_minPCT = 0.05, qc_maxMTratio = 0.1, nc_nNet = 10,
                          nc_nCells = 500, nc_nComp = 3, nc_symmetric = FALSE, nc_scaleScores = TRUE,
                          nc_q = 0.05, td_K = 3, td_maxIter = 1e3, td_maxError = 1e-5, ma_nDim = 30, dc_minDist = 1e-5){
  # Single-cell Quality Control
  X <- scQC(X, minLibSize = qc_minLibSize, removeOutlierCells = qc_removeOutlierCells, minPCT = qc_minPCT, maxMTratio = qc_maxMTratio)
  Y <- scQC(Y, minLibSize = qc_minLibSize, removeOutlierCells = qc_removeOutlierCells, minPCT = qc_minPCT, maxMTratio = qc_maxMTratio)
  
  # Counts per million (CPM) normalization
  X <- cpmNormalization(X)
  Y <- cpmNormalization(Y)

  # Comparing gene ids.
  xNames <- rownames(X)
  yNames <- rownames(Y)

  sharedGenes <- intersect(xNames, yNames)
  nGenes <- length(sharedGenes)
  
  # Filtering out non-shared genes
  X <- X[sharedGenes,]
  Y <- Y[sharedGenes,]
  
  # Construction of gene-regulatory networks based on principal component regression (pcNet) and random subsampling.
  set.seed(1)
  xList <- makeNetworks(X = X, nCells = nc_nCells, nNet = nc_nNet, nComp = nc_nComp, scaleScores = nc_scaleScores, symmetric = nc_symmetric, q = (1-nc_q))
  set.seed(1)
  yList <- makeNetworks(X = Y, nCells = nc_nCells, nNet = nc_nNet, nComp = nc_nComp, scaleScores = nc_scaleScores, symmetric = nc_symmetric, q = (1-nc_q))

  # CANDECOMP/PARAFRAC Tensor Decomposition
      # for(M in c('I','3d','4d')){
  set.seed(1)
  tensorOut <- tensorDecomposition(xList, yList, K = td_K, maxIter = td_maxIter, maxError = td_maxError)
      # Matrix::writeMM(tensorOut$X,paste0('X_',id,'_',M,'tensor.mtx'))
      # Matrix::writeMM(tensorOut$Y,paste0('Y_',id,'_',M,'tensor.mtx'))
      # writeLines(sharedGenes, paste0('genes_',id,'_',M,'tensor.mtx'))
  
  # Split of tensor output
  tX <- as.matrix(tensorOut$X)
  tY <- as.matrix(tensorOut$Y)
  
  # Non-linear manifold alignment
      # for(A in c('O','D','P')){
  set.seed(1)
  mA <- manifoldAlignment(tX , tY, d = ma_nDim)
  rownames(mA) <- c(paste0('X_', sharedGenes),paste0('y_', sharedGenes))
      # outFile <-paste0(id,'_',M,'tensor_',A,'alignment.csv')
      # write.csv(mA, outFile)
  
  # Differential coexpression testing
  dC <- dCoexpression(manifoldOutput = mA, minDist = dc_minDist)
      # write.csv(dC, paste0('dCoex_',id,'_',M,'tensor_',A,'alignment.csv'))
      # }
      # }
  
  # Return preparation
  outputResult <- list()
  outputResult$tensorNetworks <- tensorOut
  outputResult$manifoldAlignment <- mA
  outputResult$diffCoexpression <- dC

  # Return
  return(outputResult)
}
