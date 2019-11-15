#' @export makeNetworks
#' @importFrom Matrix Matrix
#' @title makeNetworks
#' @description Generates gene regulatory networks based on principal components regression for a set of subsampled cells.
#' @param X A filtered and normlized gene expression matrix with cells as columns and genes as rows.
#' @param nNet An integer value. The number of networks based on principal components regression to generate.
#' @param nCells An integer value. The number of cells to subsample each time to generate a network.
#' @param nComp An integer value. The number of principal components in PCA to generate the networks. Should be greater than 2 and lower than the total number of genes.
#' @param scaleScores A boolean value (TRUE/FALSE), if TRUE, the weights will be normalized such that the maximum absolute value is 1.
#' @param symmetric A boolean value (TRUE/FALSE), if TRUE, the weights matrix returned will be symmetric.
#' @param q A decimal value between 0 and 1. Represent the cut-off threshold of top q\% relationships to be returned.
#' @return A list of gene regulatory networks in dgCMatrix format. Each one corresponding to one group of subsampled cells.
#' @references ...
#'
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
