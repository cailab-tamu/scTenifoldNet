#' @export dCoexpression
#' @importFrom stats dist

dCoexpression <- function(manifoldOutput, nGenes, geneList, q = 0.05){
  dMetric <- sapply(seq_len(nGenes), function(G){
    X <- manifoldOutput[G,]
    Y <- manifoldOutput[(G+nGenes),]
    I <- rbind(X,Y)
    O <- dist(I)
    O <- as.numeric(O)
    return(O)
  })
  geneList[dMetric > quantile(dMetric, (1-q))]
}
