#' @export makeNetworks
#' @importFrom Matrix Matrix
#' @title Computes gene regulatory networks for subsamples of cells based on principal component regression.
#' @description This function computes \code{nNet} gene regulatory networks for a randomly selected subsample of \code{nCells} cells based on principal component regression (PCR), a technique based on principal component analysis. In PCR, the outcome variable is regressed over a \code{nComp} number of for principal components computed from a set of covariates to estimate the unknown regression coefficients in the model. \code{pcNet} function computes the PCR coefficients for each gene one at a time using all the others as covariates, to construct an all by all gene regulatory network.
#' @param X A filtered and normlized gene expression matrix with cells as columns and genes as rows.
#' @param nNet An integer value. The number of networks based on principal components regression to generate.
#' @param nCells An integer value. The number of cells to subsample each time to generate a network.
#' @param nComp An integer value. The number of principal components in PCA to generate the networks. Should be greater than 2 and lower than the total number of genes.
#' @param scaleScores A boolean value (\code{TRUE/FALSE}), if \code{TRUE}, the weights will be normalized such that the maximum absolute value is 1.
#' @param symmetric A boolean value (\code{TRUE/FALSE}), if \code{TRUE}, the weights matrix returned will be symmetric.
#' @param q A decimal value between 0 and 1. Represent the cut-off threshold of top q\% relationships to be returned.
#' @return A list with \code{nNet} gene regulatory networks in dgCMatrix format. Each one computed from a randomly selected subsample of \code{nCells} cells.
#' @references 
#' \itemize{
#' \item Gill, Ryan, Somnath Datta, and Susmita Datta. "dna: An R package for differential network analysis." Bioinformation 10.4 (2014): 233.
#' \item Jolliffe, Ian T. "A note on the use of principal components in regression." Journal of the Royal Statistical Society: Series C (Applied Statistics) 31.3 (1982): 300-303.
#' \item Massy, William F. "Principal components regression in exploratory statistical research." Journal of the American Statistical Association 60.309 (1965): 234-256.
#' }
#' @details Principal component regression may be broadly divided into three major steps: \enumerate{
#' \item Perform PCA on the observed covariates data matrix to obtain \code{nComp} number of the principal components.
#' \item Regress the observed vector of outcomes on the selected principal components as covariates, using ordinary least squares regression to get a vector of estimated regression coefficients
#' \item Transform this vector back to the scale of the actual covariates, using the eigenvectors corresponding to the selected principal components to get the final PCR estimator for estimating the regression coefficients characterizing the original model.
#' }
#' 
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
#' mnOutput <- makeNetworks(X = X,
#'                          nNet = 3, 
#'                          nCells = 500, 
#'                          nComp = 3, 
#'                          scaleScores = TRUE, 
#'                          symmetric = FALSE, 
#'                          q = 0.95
#'                          )
#' 
#' # Verifying the class
#' class(mnOutput)
#' 
#' # Verifying the number of networks
#' length(mnOutput)
#' 
#' # Veryfying the dimention of the networks
#' lapply(mnOutput,dim)
#' 
#' # Single-cell gene regulatory networks
#' mnOutput[[1]][1:10,1:10]
#' mnOutput[[2]][1:10,1:10]
#' mnOutput[[3]][1:10,1:10]

makeNetworks <- function(X, nNet = 10, nCells = 500, nComp = 3, scaleScores = TRUE, symmetric = FALSE, q = 0.95){
  geneList <- rownames(X)
  nGenes <- length(geneList)
  nCol <- ncol(X)
  if(nGenes > 0){
    pbapply::pbsapply(seq_len(nNet), function(W){
      Z <- sample(x = seq_len(nCol), size = nCells, replace = TRUE)
      Z <- as.matrix(X[,Z])
      Z <- Z[apply(Z,1,sum) > 0,]
      if(nComp > 2 & nComp < nGenes){
        Z <- pcNet(Z, nComp = nComp, scaleScores = scaleScores, symmetric = symmetric, q = q, verbose = FALSE)  
      } else {
        stop('nComp should be greater than 2 and lower than the total number of genes')
      }
      O <- matrix(data = 0, nrow = nGenes, ncol = nGenes)
      rownames(O) <- colnames(O) <- geneList
      O[rownames(Z), colnames(Z)] <- as.matrix(Z)
      O <- as(O, 'dgCMatrix')
      return(O)
    })  
  } else {
    stop('Gene names are required')
  }
}
