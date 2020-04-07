#' @export manifoldAlignment
#' @importFrom RSpectra eigs
#' @title Performs non-linear manifold alignment of two gene regulatory networks. 
#' @description Build comparable low-dimensional features for two weight-averaged denoised single-cell gene regulatory networks. Using a non-linear network embedding method \code{manifoldAlignment } aligns two gene regulatory networks and finds the structural similarities between them. This function is a wrapper of the \code{Python} code provided by Vu et al., (2012) at https://github.com/all-umass/ManifoldWarping. 
#' @param X A gene regulatory network.
#' @param Y A gene regulatory network.
#' @param d The dimension of the low-dimensional feature space.
#' @return A low-dimensional projection for two the two gene regulatory networks used as input. The output is a labeled matrix with two times the number of shared genes as rows ( X_ genes followed by Y_ genes in the same order) and \code{d} number of columns.
#' @references \itemize{
#' \item Vu, Hoa Trong, Clifton Carey, and Sridhar Mahadevan. "Manifold warping: Manifold alignment over time." Twenty-Sixth AAAI Conference on Artificial Intelligence. 2012.
#' \item Wang, Chang, and Sridhar Mahadevan. "A general framework for manifold alignment." 2009 AAAI Fall Symposium Series. 2009.
#' }
#' @details  Manifold alignment builds connections between two or more disparate data sets by aligning their underlying manifolds and provides knowledge transfer across the data sets. For further information please see: Wang et al., (2009)
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
#' # Computing 3 single-cell gene regulatory networks each one from a subsample of 500 cells
#' xNetworks <- makeNetworks(X = X,
#'                          nNet = 3, 
#'                          nCells = 500, 
#'                          nComp = 3, 
#'                          scaleScores = TRUE, 
#'                          symmetric = FALSE, 
#'                          q = 0.95
#'                          )
#' 
#' # Computing a K = 3 CANDECOMP/PARAFAC (CP) Tensor Decomposition 
#' tdOutput <- tensorDecomposition(xNetworks, K = 3, maxError = 1e5, maxIter = 1e3)
#'
#' \dontrun{
#' # Computing the alignment
#' # For this example, we are using the same input, the match should be perfect. 
#' maOutput <- manifoldAlignment(tdOutput$X, tdOutput$X)
#' 
#' # Separating the coordinates for each gene
#' X <- maOutput[grepl('X_', rownames(maOutput)),]
#' Y <- maOutput[grepl('Y_', rownames(maOutput)),]
#' 
#' # Plotting
#' # X Points
#' plot(X, pch = 16)
#' 
#' # Y Points
#' points(Y, col = 'red')
#' 
#' # Legend
#' legend('topright', legend = c('X', 'Y'), 
#'        col = c('black', 'red'), bty = 'n', 
#'        pch = c(16,1), cex = 0.7)
#' }

manifoldAlignment <- function(X, Y, d = 30){
  sharedGenes <- intersect(rownames(X), rownames(Y))
  X <- X[sharedGenes, sharedGenes]
  Y <- Y[sharedGenes, sharedGenes]
  L <- diag(length(sharedGenes))
  wX <- X+1
  wY <- Y+1
  wXY <- 0.9 * (sum(wX) + sum(wY)) / (2 * sum(L)) * L
  W <- rbind(cbind(wX, wXY), cbind(t(wXY), wY))
  W <- -W
  diag(W) <- 0
  diag(W) <- -apply(W, 2, sum)
  E <- suppressWarnings(RSpectra::eigs(W, d*2, 'SR'))
  E$values <- suppressWarnings(as.numeric(E$values))
  E$vectors <- suppressWarnings(apply(E$vectors,2,as.numeric))
  newOrder <- order(E$values)
  E$values <- E$values[newOrder]
  E$vectors <- E$vectors[,newOrder]
  E$vectors <- E$vectors[,E$values > 1e-8]
  alignedNet <- E$vectors[,1:30]
  colnames(alignedNet) <- paste0('NLMA ', seq_len(d))
  rownames(alignedNet) <- c(paste0('X_', sharedGenes), paste0('Y_', sharedGenes))
  return(alignedNet)
}
