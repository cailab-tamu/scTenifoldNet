library(Matrix)
library(scTenifoldNet)

setwd('/data/dcosorioh/DF')

X <- readMM('unstimulated/matrix.mtx')
rownames(X) <- read.csv('unstimulated/genes.tsv', sep = '\t', header = FALSE, stringsAsFactors = FALSE)[,2]
colnames(X) <- readLines('unstimulated/barcodes.tsv')
X <- X[rowMeans(X != 0) > 0.05,]

Y <- readMM('stimulated/matrix.mtx')
rownames(Y) <- read.csv('stimulated/genes.tsv', sep = '\t', header = FALSE, stringsAsFactors = FALSE)[,2]
colnames(Y) <- readLines('stimulated/barcodes.tsv')
Y <- Y[rowMeans(Y != 0) > 0.05,]

#DF <- scTenifoldNet(X, Y)
#save(DF, file = 'DF_1.RData')

X <- X[!grepl('^RPL|^RPS|MT-',rownames(X)),]
Y <- Y[!grepl('^RPL|^RPS|MT-',rownames(Y)),]
 
DF <- scTenifoldNet(X, Y)
save(DF, file = 'DF_2.RData')