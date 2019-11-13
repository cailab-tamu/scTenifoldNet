#' @export scTenifoldNet
#' @title scTenifoldNet
#' @importFrom methods as
#' @importFrom rTensor as.tensor cp
#' @description ...
#' @param X ...
#' @param Y ...
#' @param qc_minLibSize ...
#' @param qc_removeOutlierCells ...
#' @param qc_minPCT ...
#' @param qc_maxMTratio ...
#' @param nNet ...
#' @param nCells ...
#' @param nComp ...
#' @param symmetric ...
#' @param scaleScores ...
#' @param q ...
#' @return ...
#' @references ...

scTenifoldNet <- function(X, Y, qc_minLibSize = 1000, qc_removeOutlierCells = TRUE,
                          qc_minPCT = 0.05, qc_maxMTratio = 0.1, nNet = 10, 
                          nCells = 500, nComp = 3, symmetric = FALSE, scaleScores = TRUE, 
                          q = 0.05){
  
  X <- scQC(X,minLibSize = qc_minLibSize, removeOutlierCells = qc_removeOutlierCells, minPCT = qc_minPCT, maxMTratio = qc_maxMTratio)
  X <- cpmNormalization(X)
  
  Y <- scQC(Y)
  Y <- cpmNormalization(Y)
  
  xNames <- rownames(X)
  yNames <- rownames(Y)
  
  sharedGenes <- intersect(xNames, yNames)
  nGenes <- length(sharedGenes)
  
  X <- X[sharedGenes,]
  Y <- Y[sharedGenes,]
  
  set.seed(1)
  xList <- makeNetworks(X = X, nCells = nCells, nNet = nNet, nComp = nComp, scaleScores = scaleScores, symmetric = symmetric, q = (1-q))
  set.seed(1)
  yList <- makeNetworks(X = Y, nCells = nCells, nNet = nNet, nComp = nComp, scaleScores = scaleScores, symmetric = symmetric, q = (1-q))
  
  # for(M in c('I','3d','4d')){
  set.seed(1)
  tensorOut <- tensorDecomposition(xList, yList, d = 3)
  # Matrix::writeMM(tensorOut$X,paste0('X_',id,'_',M,'tensor.mtx'))
  # Matrix::writeMM(tensorOut$Y,paste0('Y_',id,'_',M,'tensor.mtx'))
  # writeLines(sharedGenes, paste0('genes_',id,'_',M,'tensor.mtx'))
  tX <- as.matrix(tensorOut$X)
  tY <- as.matrix(tensorOut$Y)
  # for(A in c('O','D','P')){
  set.seed(1)
  mA <- manifoldAlignment(tX , tY, d = 30)
  rownames(mA) <- c(paste0('X_', sharedGenes),paste0('y_', sharedGenes))
  # outFile <-paste0(id,'_',M,'tensor_',A,'alignment.csv')
  # write.csv(mA, outFile)
  dC <- dCoexpression(mA, length(sharedGenes), sharedGenes)
  # write.csv(dC, paste0('dCoex_',id,'_',M,'tensor_',A,'alignment.csv'))
  # }
  # }
  outputResult <- list()
  outputResult$Networks <- tensorOut
  outputResult$ManifoldAlignment <- mA
  outputResult$diffCoexpression <- dC
  
  return(outputResult)
}