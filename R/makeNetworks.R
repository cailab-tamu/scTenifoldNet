#' @export makeNetworks
#' @importFrom Matrix Matrix
#' @title makeNetworks
#' @description Generate principal components gene regulatory networks for a list of subsampled cells.
#' @param X Normlized gene expression matrix with cells as columns and genes (symbols) as rows.
#' @param nNet The number of principal components regression based networks to generate.
#' @param nCells The number of cells to subsample each time.
#' @param nComp The number of principal components in PCA.
#' @param scaleScores A boolean value (TRUE/FALSE), if TRUE, the weights will be normalized such that the maximum absolute value is 1.
#' @param symmetric A boolean value (TRUE/FALSE), if TRUE, the weights matrix returned will be symmetric.
#' @param q The threshold that only remaining the relationships with the top q% absolute value in the weights matrix.
#' @return A list of gene regulatory networks. Each one corresponding to one group of subsampled cells.
#' @references ...
#'
makeNetworks <- function(X, nNet = 10, nCells = 500, nComp = 3, scaleScores = TRUE, symmetric = FALSE, q = 0.95){
  sapply(seq_len(nNet), function(W){
    Z <- sample(x = seq_len(ncol(X)), size = nCells, replace = TRUE)
    Z <- as.matrix(X[,Z])
    Z <- Z[apply(Z,1,sum) > 0,]
    Z <- pcNet(Z, nCom = nComp, scaleScores = scaleScores, symmetric = symmetric, q = q)
    geneList <- rownames(X)
    nGenes <-length(geneList)
    O <- matrix(data = 0, nrow = nGenes, ncol = nGenes)
    rownames(O) <- colnames(O) <- geneList
    O[rownames(Z), colnames(Z)] <- as.matrix(Z)
    O <- Matrix::Matrix(O)
    return(O)
  })
}
