#' @export makeNetworks
#' @importFrom Matrix Matrix
makeNetworks <- function(X, nNet = 10, nCells = 1000, nComp = 3, scaleScores = TRUE, symmetric = FALSE, q = 0.95){
  sapply(seq_len(nNet), function(W){
    Z <- sample(x = colnames(X), size = nCells, replace = TRUE)
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
