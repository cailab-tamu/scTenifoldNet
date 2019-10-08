library(PCrTdMa)
library(Matrix)

method <- 'PCR'
nRun <- 10
nCells <- c(10, 50, 100, 500, 1000, 2000, 3000)
nRun <- seq_len(nRun)
nCells <- nCells

experimentDesign <- expand.grid(method, nRun, nCells)
experimentDesign <- as.data.frame.array(experimentDesign)
dir.create(paste0('networks/',method))

inputData <- readMM('data/CTL.mtx')
rownames(inputData) <- readLines('data/geneList.txt')

set.seed(1)
apply(experimentDesign,1,function(X){
  newFile <- paste0('networks/', method, '/', X[1], '_run', as.character(X[2]), '_',as.character(X[3]), 'cells', collapse = '')
  sCells <- sample(seq_len(ncol(inputData)), size = as.numeric(X[3]))
  tData <- inputData[,sCells]
  gList <- rownames(tData)
  nGenes <- length(gList)
  outM <- matrix(0, nrow = nGenes, ncol = nGenes)
  rownames(outM) <- colnames(outM) <- gList
  tData <- tData[apply(tData!=0,1,sum) > 0,]
  tData <- as.matrix(tData)
  tData <- pcNet(tData)
  tData <- round(tData,2)
  outM[rownames(tData),colnames(tData)] <- as.matrix(tData)
  outM <- as(outM, 'dgCMatrix')
  writeMM(outM, newFile)
})

