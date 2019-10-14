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

dCoexpression <- function(manifoldOutput, nGenes, geneList, q = 0.05){
  dMetric <- sapply(seq_len(nGenes), function(G){
    X <- manifoldOutput[G,]
    Y <- manifoldOutput[(G+nGenes),]
    I <- rbind(X,Y)
    O <- dist(I)
    O <- as.numeric(O)
    return(O)
  })
  
  pValues <- pchisq(abs(scale(dMetric)), df = 1, lower.tail = FALSE)
  pAdjusted <- p.adjust(pValues, method = 'fdr')
  dOut <- data.frame(
    gene = geneList, 
    distance = dMetric,
    p.value = pValues,
    p.adj = pAdjusted
  )
  return(dOut)
}
