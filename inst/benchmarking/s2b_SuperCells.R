library(Matrix)
library(SuperCell)
library(pbapply)
library(scTenifoldNet)

method <- 'SUPERCELL_PC'
nRun <- 10
nCells <- c(10, 50, 100, 500, 1000, 2000, 3000)
nRun <- seq_len(nRun)
nCells <- nCells

experimentDesign <- expand.grid(method, nRun, nCells)
experimentDesign <- as.data.frame.array(experimentDesign)
dir.create(paste0('networks/',method))

inputData <- readMM('data/CTL.mtx')
rownames(inputData) <- readLines('data/geneList.txt')

O <- SuperCell::SCimplify(inputData)
O <- SuperCell::supercell_GE(inputData, O$membership)
#O <- log1p(t(t(O)/colSums(O)) * 1e4)
O <- preprocessCore::normalize.quantiles(as.matrix(O))

dimnames(O) <- dimnames(inputData)
inputData <- O


set.seed(1)
apply(experimentDesign,1,function(X){
  newFile <- paste0('networks/', method, '/', X[1], '_run', as.character(X[2]), '_',as.character(X[3]), 'cells', collapse = '')
  newFile <- gsub('[[:space:]]+','',newFile)
  sCells <- sample(seq_len(ncol(inputData)), size = as.numeric(X[3]), replace = TRUE)
  tData <- inputData[,sCells]
  gList <- rownames(tData)
  nGenes <- length(gList)
  outM <- matrix(0, nrow = nGenes, ncol = nGenes)
  rownames(outM) <- colnames(outM) <- gList
  tData <- tData[apply(tData!=0,1,sum) > 0,]
  tData <- as.matrix(tData)
  tData <- cor(t(tData))
  #tData <- pcNet(tData)
  tData <- round(tData,2)
  outM[rownames(tData),colnames(tData)] <- as.matrix(tData)
  outM <- as(outM, 'dgCMatrix')
  writeMM(outM, newFile)
})


q <- seq(0,1,0.05)
metricOutput <- pblapply(q, function(Q){
  Acurracy <- lapply(nCells, function(X){
    fileList <- list.files('networks/SUPERCELL_PC/', full.names = TRUE)
    fileList <- fileList[grepl(paste0(X,'cells'), fileList)]
    fileContent <- lapply(fileList, readMM)
    fileContent <- lapply(fileContent, function(X){
      X <- as.matrix(X)
      X <- X/max(abs(X))
      X <- (X - 0)
      X[abs(X) < quantile(abs(X), Q)] <-  0
      diag(X) <- 1
      dataMatrix <- X
      TP <- sum(dataMatrix[1:40,1:40] > 0) + sum(dataMatrix[41:98,41:98] > 0)
      FP <- sum(dataMatrix[1:40,41:98] > 0) + sum(dataMatrix[41:98, 1:40] > 0)
      TN <- sum(dataMatrix[1:40,41:98] < 0) + sum(dataMatrix[41:98, 1:40] < 0)
      FN <- sum(dataMatrix[1:40,1:40] < 0) + sum(dataMatrix[41:98,41:98] < 0)
      ACC <- round((TP+TN)/(TP+TN+FN+FP),2)
      return(ACC)
    })
    unlist(fileContent)
  })
  ReCall <- lapply(nCells, function(X){
    fileList <- list.files('networks/SUPERCELL_PC/', full.names = TRUE)
    fileList <- fileList[grepl(paste0(X,'cells'), fileList)]
    fileContent <- lapply(fileList, readMM)
    fileContent <- lapply(fileContent, function(X){
      X <- as.matrix(X)
      X <- X/max(abs(X))
      X <- (X - 0)
      X[abs(X) < quantile(abs(X), Q)] <-  0
      diag(X) <- 1
      dataMatrix <- X
      TP <- sum(dataMatrix[1:40,1:40] > 0) + sum(dataMatrix[41:98,41:98] > 0)
      FP <- sum(dataMatrix[1:40,41:98] > 0) + sum(dataMatrix[41:98, 1:40] > 0)
      TN <- sum(dataMatrix[1:40,41:98] < 0) + sum(dataMatrix[41:98, 1:40] < 0)
      FN <- sum(dataMatrix[1:40,1:40] < 0) + sum(dataMatrix[41:98,41:98] < 0)
      REC <- round((TP)/((40*40)+(58*58)),2)
      return(REC)
    })
    unlist(fileContent)
  })
  
  accMean <- unlist(lapply(Acurracy, mean))
  accSD <- unlist(lapply(Acurracy, sd))
  recallMean <- unlist(lapply(ReCall, mean))
  recallSD <- unlist(lapply(ReCall, sd))
  
  outputMetric <- NULL
  outputMetric$q <- rep(Q, length(nCells))
  outputMetric$nCells <- nCells
  outputMetric$accLB <- accMean-accSD
  outputMetric$acc <- accMean
  outputMetric$accUB <- accMean+accSD
  outputMetric$recallLB <- recallMean - recallSD
  outputMetric$recall <- recallMean
  outputMetric$recallUB <- recallMean + recallSD
  
  outputMetric <- as.data.frame(outputMetric)
  outputMetric
})

metricOutput <- do.call(rbind.data.frame, metricOutput)
metricOutput <- round(metricOutput,3)
write.csv(metricOutput, row.names = FALSE, file = 'metrics/SUPERCELL_PC.csv')


