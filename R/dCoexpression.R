#' @export dCoexpression
#' @importFrom stats dist pnorm p.adjust
#' @importFrom MASS boxcox
#' @title Evaluates gene differential coexpression based on manifold alignment distances. 
#' @description Using the output of the non-linear manifold alignment, this function computes the Euclidean distance between the coordinates for the same gene in both conditions. Calculated distances are then transformed using Box-Cox power transformation, and standardized to ensure normality. P-values are assigned following the standard normal distribution. 
#' @param manifoldOutput A matrix. The output of the non-linear manifold alignment,  a labeled matrix with two times the number of shared genes as rows (X_ genes followed by Y_ genes in the same order) and \code{d} number of columns.
#' @param minDist A decimal value. Defines the cut-off threshold of distance to limit the testing to genes that show, at least \code{minDist} deviation.
#' @return A data frame with 5 columns as follows: \itemize{
#' \item \code{gene} A character vector with the gene id identified from the \code{manifoldAlignment} output.
#' \item \code{distance} A numeric vector of the Euclidean distance computed between the coordinates of the same gene in both conditions.
#' \item \code{Z} A numeric vector of the Z-scores computed after Box-Cox power transformation.
#' \item \code{p.value} A numeric vector of the p-values associated to the Z-scores, probabilities are \eqn{P[X > x]}
#' \item \code{p.adj} A numeric vector of corrected p-values using Benjamini & Hochberg (1995) FDR.
#' }
#' @references \itemize{
#' \item Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics, 29, 1165-1188. doi: 10.1214/aos/1013699998.
#' }
#' 
dCoexpression <- function(manifoldOutput, minDist = 1e-5){
  
  geneList <- rownames(manifoldOutput)
  geneList <- geneList[grepl('^X_', geneList)]
  geneList <- gsub('^X_','', geneList)
  nGenes <- length(geneList)
  
  eGenes <- nrow(manifoldOutput)/2
  
  eGeneList <- rownames(manifoldOutput)
  eGeneList <- eGeneList[grepl('^Y_', eGeneList)]
  eGeneList <- gsub('^Y_','', eGeneList)
  
  if(nGenes != eGenes){
    stop('Number of identified and expected genes are not the same')
  }
  if(!all(eGeneList == geneList)){
    stop('Genes are not ordered as expected. X_ genes should be followed by Y_ genes in the same order')
  }
  
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
  
  BC <- boxcox(dMetric~1, lambda = seq(-3,3,1/1000), plotit = FALSE)
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
  dOut <- as.data.frame.array(dOut)
  return(dOut)
}

