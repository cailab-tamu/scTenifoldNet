library(scTenifoldNet)
library(Matrix)

X <- read.csv(file = 'X.txt', header = FALSE, sep = '\t')
X <- as.matrix(X)
rownames(X) <- readLines('genelist.txt')


Y <- read.csv(file = 'Y.txt', header = FALSE, sep = '\t')
Y <- as.matrix(Y)
rownames(Y) <- readLines('genelist.txt')

output <- scTenifoldNet(X = X, Y = Y,
                        nc_nNet = 10, nc_nCells = 500,
                        td_K = 3, qc_minLibSize = 30,
                        dc_minFC = 0)

output$diffRegulation[1:5,]