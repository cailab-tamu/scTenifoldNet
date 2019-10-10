library(Matrix)
library(rTensor)
library(PCrTdMa)
source('../../R/avgNetworks.R')

X <- readMM('data/CTL.mtx')
rownames(X) <- readLines('data/geneList.txt')

nCells <- c(500, 1000, 2000, 3000)
q <- c(0.95, 0.9, 0.85, 0.8)
nNets <- c(5, 10, 15, 20, 50, 100)

experimentDesign <- expand.grid(nCells, nNets, q)
experimentDesign <- experimentDesign[order(experimentDesign$Var1),]

dir.create('networks/TENSOR')

apply(experimentDesign,1,function(P){
  avgFile <- paste0('networks/TENSOR/avg_',P[1],'cells_', P[2], 'nets_', P[3],'q.mtx')
  tsrFile <- paste0('networks/TENSOR/tsr_',P[1],'cells_', P[2], 'nets_', P[3],'q.mtx')

  nCells <- as.numeric(P[1])
  nNets <- as.numeric(P[2])
  q <- as.numeric(P[3])

  set.seed(1)
  netX <- makeNetworks(X, nNet = nNets, nCells = nCells, q = q)
  avgX <- round(avgNetworks(netX),3)
  writeMM(avgX, avgFile)

  Tensor <- array(dim = c(100,100,1,nNets))
  for(i in seq_along(netX)){
    Tensor[,,,i] <- as.matrix(netX[[i]])
  }
  Tensor <- as.tensor(Tensor)
  Tensor <- cp(Tensor, num_components = 3, max_iter = 1000)

  tsrX <- NULL
  for(i in seq_along(netX)){
    tsrX[[i]] <- Tensor$est@data[,,,i]
  }
  tsrX <- round(avgNetworks(tsrX),3)
  tsrX <- as(tsrX, 'dgCMatrix')
  writeMM(tsrX, tsrFile)
})
