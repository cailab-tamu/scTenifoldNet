#' @export dRegulation
#' @importFrom stats dist pchisq p.adjust qqnorm
#' @importFrom MASS boxcox
#' @title Evaluates gene differential regulation based on manifold alignment distances. 
#' @description Using the output of the non-linear manifold alignment, this function computes the Euclidean distance between the coordinates for the same gene in both conditions. Calculated distances are then transformed using Box-Cox power transformation, and standardized to ensure normality. P-values are assigned following the chi-square distribution over the fold-change of the squared distance computed with respect to the expectation. 
#' @param manifoldOutput A matrix. The output of the non-linear manifold alignment,  a labeled matrix with two times the number of shared genes as rows (X_ genes followed by Y_ genes in the same order) and \code{d} number of columns.
#' @return A data frame with 6 columns as follows: \itemize{
#' \item \code{gene} A character vector with the gene id identified from the \code{manifoldAlignment} output.
#' \item \code{distance} A numeric vector of the Euclidean distance computed between the coordinates of the same gene in both conditions.
#' \item \code{Z} A numeric vector of the Z-scores computed after Box-Cox power transformation.
#' \item \code{FC} A numeric vector of the FC computed with respect to the expectation.
#' \item \code{p.value} A numeric vector of the p-values associated to the fold-changes, probabilities are asigned as \eqn{P[X > x]} using the Chi-square distribution with one degree of freedom.
#' \item \code{p.adj} A numeric vector of adjusted p-values using Benjamini & Hochberg (1995) FDR correction.
#' }
#' @references \itemize{
#' \item Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics, 29, 1165-1188. doi: 10.1214/aos/1013699998.
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
#' xNetworks <- makeNetworks(X = qcOutput,
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
#' # Evaluating the difference in regulation
#' dcOutput <- dRegulation(maOutput, minFC = 0)
#' head(dcOutput)
#' 
#' # Plotting
#' # If FDR < 0.05, the gene will be colored in red.
#' geneColor <- ifelse(dcOutput$p.adj < 0.05, 'red', 'black')
#' qqnorm(dcOutput$Z, main = 'Standardized Distance', pch = 16, col = geneColor)
#' qqline(dcOutput$Z)
#' }

dRegulation <- function(manifoldOutput){
  
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
  
  ### BOX-COX
  lambdaValues <- seq(-2,2,length.out = 1000)
  lambdaValues <- lambdaValues[lambdaValues != 0]
  BC <- try(MASS::boxcox(dMetric~1, plot=FALSE, lambda = lambdaValues), silent = TRUE)
  if(class(BC) == 'try-error'){
    nD <- dMetric
  } else {
    BC <- BC$x[which.max(BC$y)]
    if(BC < 0){
      nD <- 1/(dMetric ^ BC)
    } else {
      nD <- dMetric ^ BC
    }
  }
  Z <- scale(nD)
  E <- mean(dMetric^2)
  FC <- dMetric^2/E
  pValues <- pchisq(q = FC,df = 1,lower.tail = FALSE)
  pAdjusted <- p.adjust(pValues, method = 'fdr')
  dOut <- data.frame(
    gene = geneList, 
    distance = dMetric,
    Z = Z,
    FC = FC,
    p.value = pValues,
    p.adj = pAdjusted
  )
  dOut <- dOut[order(dOut$p.value),]
  dOut <- as.data.frame.array(dOut)
  return(dOut)
}

