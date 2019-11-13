#' @export dCoexpression
#' @importFrom stats dist
#' @title dCoexpression
#' @description ...
#' @param manifoldOutput ...
#' @param nGenes ...
#' @param geneList ...
#' @param q ...
#' @return ...
#' @references ...

dCoexpression <- function(manifoldOutput, nGenes, geneList){
  dMetric <- sapply(seq_len(nGenes), function(G){
    X <- manifoldOutput[G,]
    Y <- manifoldOutput[(G+nGenes),]
    I <- rbind(X,Y)
    O <- dist(I)
    O <- as.numeric(O)
    return(O)
  })
  
  dMetric[is.na(dMetric)] <- 0
  sDist <- scale(dMetric)
  
  pValues <- pnorm(sDist, lower.tail = FALSE)
  pAdjusted <- p.adjust(pValues, method = 'fdr')
  dOut <- data.frame(
    gene = geneList, 
    distance = dMetric,
    Z = sDist,
    p.value = pValues,
    p.adj = pAdjusted
  )
  dOut <- dOut[order(dOut$p.value),]
  return(dOut)
}
