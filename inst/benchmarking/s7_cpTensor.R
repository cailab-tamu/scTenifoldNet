library(Matrix)
library(rTensor)
source('../../R/avgNetworks.R')

X <- readMM('data/CTL.mtx')
rownames(X) <- readLines('data/geneList.txt')

nCells <- c(500, 1000, 2000, 3000)
q <- c(0.95, 0.9, 0.85, 0.8)
nNets <- c(5, 10, 15, 20, 50, 100)

experimentDesign <- expand.grid(nCells, nNets, q)

dir.create('networks/TENSOR')

# apply(experimentDesign,1,function(P){
#   avgFile <- paste0('networks/TENSOR/avg_',P[1],'cells_', P[2], 'nets_', P[3],'q.mtx')
#   tsrFile <- paste0('networks/TENSOR/tsr_',P[1],'cells_', P[2], 'nets_', P[3],'q.mtx')
#
#   nCells <- as.numeric(P[1])
#   nNets <- as.numeric(P[2])
#   q <- as.numeric(P[3])
#
#   netX <- makeNetworks(X, nNet = nNets, nCells = nCells, q = q)
#   avgX <- round(avgNetworks(netX),3)
#   writeMM(avgX, avgFile)
#
#   Tensor <- array(dim = c(100,100,1,nNets))
#   for(i in seq_along(netX)){
#     Tensor[,,,i] <- as.matrix(netX[[i]])
#   }
#   Tensor <- as.tensor(Tensor)
#   Tensor <- cp(Tensor, num_components = 3, max_iter = 1000)
#
#   tsrX <- NULL
#   for(i in seq_along(netX)){
#     tsrX[[i]] <- Tensor$est@data[,,,i]
#   }
#   tsrX <- round(avgNetworks(tsrX),3)
#   tsrX <- as(tsrX, 'dgCMatrix')
#   writeMM(tsrX, tsrFile)
# })

Accuracy <- function(X){
  c1p <- sum(X[1:40,1:40] > 0)
  c2p <- sum(X[41:98,41:98] > 0)
  c1n <- sum(X[1:40,41:98] < 0)
  c2n <- sum(X[41:98,1:40] < 0)
  (c1p+c1n+c2p+c2n)/sum(X[1:98,1:98] != 0)
}

Recall <- function(X){
  c1p <- sum(X[1:40,1:40] > 0)
  c2p <- sum(X[41:98,41:98] > 0)
  (c1p+c2p)/((40*40)+(58*58))
}

PR <- t(apply(experimentDesign,1,function(P){
  avgFile <- paste0('networks/TENSOR/avg_',P[1],'cells_', P[2], 'nets_', P[3],'q.mtx')
  tsrFile <- paste0('networks/TENSOR/tsr_',P[1],'cells_', P[2], 'nets_', P[3],'q.mtx')

  avgNet <- readMM(avgFile)
  avgA <- Accuracy(avgNet)
  avgR <- Recall(avgNet)
  tsrNet <- readMM(tsrFile)
  tsrA <- Accuracy(tsrNet)
  tsrR <- Recall(tsrNet)
  c(avgA, avgR, tsrA, tsrR)
}))

PR <- cbind(experimentDesign, PR)
colnames(PR) <- c('nCells', 'nNets', 'q', 'avgA', 'avgR', 'tsrA', 'tsrR')


par(mfrow=c(1,2))
plot(0, xlim = c(0,100), ylim = c(0,1), type = 'n', ylab = 'Accuracy')
abline(h=0.5, lty=3, col='gray90')
for(i in c(0.95,0.9)){
  tPR <- PR[PR$nCells == 500 & PR$q == i,]
  points(tPR$nNets, tPR$tsrA, ylim=c(0,1), type='b', lty = 1)
  points(tPR$nNets, tPR$avgA, ylim=c(0,1), type='b', lty = 2)
}

plot(0, xlim = c(0,100), ylim = c(0,1), type = 'n', ylab = 'Recall')
abline(h=0.5, lty=3, col='gray90')
for(i in c(500,1000,2000,3000)){
  tPR <- PR[PR$nCells == 500 & PR$q == 0.9,]
  points(tPR$nNets, tPR$tsrR, ylim=c(0,1), type='b', lty = 1)
  points(tPR$nNets, tPR$avgR, ylim=c(0,1), type='b', lty = 2)
}
