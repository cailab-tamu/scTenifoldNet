#' @export dCoexpression
#' @importFrom stats dist pnorm p.adjust
#' @title Evaluates gene differential coexpression based on manifold alignment distances. 
#' @description Using the output of the non-linear manifold alignment, this function computes the Euclidean distance between the coordinates for the same gene in both conditions. Calculated distances are then transformed using Box-Cox power transformation, and standardized to ensure normality. P-values are assigned following the standard normal distribution. 
#' @param manifoldOutput ...
#' @param nGenes ...
#' @param geneList ...
#' @param minDist ...
#' @return ...
#' @references ...

dCoexpression <- function(manifoldOutput, nGenes, geneList, minDist = 1e-5){
  dMetric <- sapply(seq_len(nGenes), function(G){
    X <- manifoldOutput[G,]
    Y <- manifoldOutput[(G+nGenes),]
    I <- rbind(X,Y)
    O <- dist(I)
    O <- as.numeric(O)
    return(O)
  })
  
  dMetric[is.na(dMetric)] <- 0
  dMetric[dMetric < minDist] <- 0
  
  geneList <- geneList[dMetric != 0]
  dMetric <- dMetric[dMetric != 0]
  
  BC <- MASS::boxcox(dMetric~1, lambda = seq(-3,3,1/1000))
  BC <- BC$x[which.max(BC$y)]
  
  if(BC < 0){
    sDist <- scale(1/(dMetric ^ BC))
  } else {
    sDist <- scale(dMetric ^ BC)
  }
  
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
