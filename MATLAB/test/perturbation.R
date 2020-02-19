X0 <- read.csv('X0.csv', header = FALSE)
X1 <- read.csv('X1.csv', header = FALSE)

rownames(X0) <- paste0('g', 1:nrow(X0))
rownames(X1) <- paste0('g', 1:nrow(X1))

library(scTenifoldNet)
X <- as.matrix(X0)
Y <- as.matrix(X1)

outputResult <- scTenifoldNet(X = X, Y = Y, dc_minFC = 0, 
                              qc_minLibSize = 0, qc_removeOutlierCells = FALSE, qc_minPCT = 0,
                              qc_maxMTratio = 1)
outputResult$diffRegulation
