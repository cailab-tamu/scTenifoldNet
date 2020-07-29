library('Matrix')
library(igraph)
library(ggplot2)
library(patchwork)
library(RSpectra)

eX <- readMM('datasets/morphineNeurons/mock/control.mtx')
rownames(eX) <- readLines('datasets/morphineNeurons/mock/controlGenes.tsv')
eX <- eX[rowMeans(eX != 0) > 0.05,]

eY <- readMM('datasets/morphineNeurons/morphine/morphine.mtx')
rownames(eY) <- readLines('datasets/morphineNeurons/morphine/morphineGenes.tsv')
eY <- eY[rowMeans(eY != 0) > 0.05,]

sGenes <- intersect(rownames(eX), rownames(eY))

eX <- eX[sGenes,]
eY <- eY[sGenes,]

eX <- svds(as.matrix(t(eX)), 30)
eY <- svds(as.matrix(t(eY)), 30)

A <- eX$v
B <- eY$v

dE <- sapply(seq_along(sGenes), function(X){
  dist(rbind(A[X,],B[X,]))
})
names(dE) <- sGenes


O <- read.csv('results/sym10X500morphineNeuron_Itensor_Dalignment.csv')
Z <- O$distance
names(Z) <- O$gene

tGenes <- intersect(names(dE), names(Z))

DF <- data.frame(D=log1p(dE[tGenes]),MA=log10(Z[tGenes]))

A <- ggplot(DF, mapping = aes(D, MA)) +
  geom_point(col = densCols(DF$D, DF$MA, colramp = hcl.colors)) +
  theme_bw() +
  xlab(expression(log(Expression~Distance + 1))) +
  ylab(expression(log[10](Manifold~Distance))) +
  labs(title = 'Expression Distance vs. Manifold Distance', subtitle = parse(text = paste0('rho == ', round(cor(DF$D, DF$MA, method = 'sp'),3)))) +
  theme(plot.title = element_text(face = 2))

png('DComparison.png', width = 1300, height = 1300, res = 300)
A
dev.off()



eX <- readMM('results/tensorOutput/X_10X500morphineNeuron_Itensor.mtx')
rownames(eX) <- readLines('results/tensorOutput/genes_10X500morphineNeuron_Itensor.mtx')

eY <- readMM('results/tensorOutput/Y_10X500morphineNeuron_Itensor.mtx')
rownames(eY) <- readLines('results/tensorOutput/genes_10X500morphineNeuron_Itensor.mtx')

sGenes <- intersect(rownames(eX), rownames(eY))

eX <- eX[sGenes,]
eY <- eY[sGenes,]

eX <- svds(as.matrix(t(eX)), 30)
eY <- svds(as.matrix(t(eY)), 30)

A <- eX$v
B <- eY$v

dE <- sapply(seq_along(sGenes), function(X){
  dist(rbind(A[X,],B[X,]))
})
names(dE) <- sGenes


O <- read.csv('results/sym10X500morphineNeuron_Itensor_Dalignment.csv')
Z <- O$distance
names(Z) <- O$gene

tGenes <- intersect(names(dE), names(Z))

DF <- data.frame(D=log10(dE[tGenes]),MA=log10(Z[tGenes]))

A <- ggplot(DF, mapping = aes(D, MA)) +
  geom_point(col = densCols(DF$D, DF$MA, colramp = hcl.colors)) +
  theme_bw() +
  xlab(expression(log[10](Tensor~Distance))) +
  ylab(expression(log[10](Manifold~Distance))) +
  labs(title = 'Tensor Distance vs. Manifold Distance', subtitle = parse(text = paste0('rho == ', round(cor(DF$D, DF$MA, method = 'sp'),3)))) +
  theme(plot.title = element_text(face = 2))

png('NComparison.png', width = 1300, height = 1300, res = 300)
A
dev.off()
