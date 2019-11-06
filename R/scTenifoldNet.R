#' @export scTenifoldNet
#' @title scTenifoldNet
#' @importFrom methods as
#' @importFrom rTensor as.tensor cp
#' @description ...
#' @param X ...
#' @param Y ...
#' @param nNet ...
#' @param nCells ...
#' @param nComp ...
#' @param symmetric ...
#' @param scaleScores ...
#' @param q ...
#' @return ...
#' @references ...
#devtools::install_github('cailab-tamu/scTenifoldNet')
# library(scTenifoldNet)
# library(Matrix)
# 
# X <- readMM('../manuscript/datasets/neurons/mock/control.mtx')
# rownames(X) <- readLines('../manuscript/datasets/neurons/mock/controlGenes.tsv')
# X <- X[!grepl('Rik$',rownames(X)),]
# X <- X[1:500,]
# Y <- readMM('../manuscript/datasets/neurons/morphine/morphine.mtx')
# rownames(Y) <- readLines('../manuscript/datasets/neurons/morphine/morphineGenes.tsv')
# Y <- Y[!grepl('Rik$',rownames(Y)),]
# Y <- Y[1:500,]

scTenifoldNet <- function(X, Y, id, nNet = 10, nCells = 500, nComp = 3, symmetric = FALSE, scaleScores = TRUE, q = 0.05){
  
  X <- scTenifoldNet:::scQC(X)
  X <- scTenifoldNet:::CPM(X)
  X <- X[apply(X!=0,1,mean) > 0.05,]
  
  Y <- scTenifoldNet:::scQC(Y)
  Y <- scTenifoldNet:::CPM(Y)
  Y <- Y[apply(Y!=0,1,mean) > 0.05,]
  
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
  
  for(M in c('I','3d','4d')){
    set.seed(1)
    tensorOut <- tensorDecomposition(xList, yList, d = 3, type = M)
    Matrix::writeMM(tensorOut$X,paste0('X_',id,'_',M,'tensor.mtx'))
    Matrix::writeMM(tensorOut$Y,paste0('Y_',id,'_',M,'tensor.mtx'))
    writeLines(sharedGenes, paste0('genes_',id,'_',M,'tensor.mtx'))
    tX <- as.matrix(tensorOut$X)
    tY <- as.matrix(tensorOut$Y)
    for(A in c('O','D','P')){
      set.seed(1)
      mA <- manifoldAlignment(tX , tY, d = 100, type = A)
      rownames(mA) <- c(paste0('X_', sharedGenes),paste0('y_', sharedGenes))
      outFile <-paste0(id,'_',M,'tensor_',A,'alignment.csv')
      write.csv(mA, outFile)
      dC <- dCoexpression(mA, length(sharedGenes), sharedGenes)
      write.csv(dC, paste0('dCoex_',id,'_',M,'tensor_',A,'alignment.csv'))
    }
  }
}

#scTenifoldNet(X = X, Y = Y, id = '500morphineNeuron', nNet = 5, nCells = 1000)
