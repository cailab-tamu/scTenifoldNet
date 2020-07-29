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
BC <- MASS::boxcox(dE~1)
dE <- dE ^ abs(BC$x[which.max(BC$y)])
zDE <- as.vector(scale(dE))
names(zDE) <- names(dE)

O <- read.csv('results/sym10X500morphineNeuron_Itensor_Dalignment.csv')
Z <- O$Z
names(Z) <- O$gene

tGenes <- intersect(names(zDE), names(Z))

DF <- data.frame(D=(zDE[tGenes]),MA=(Z[tGenes]))

A <- ggplot(DF, mapping = aes(D, MA)) +
  geom_point(col = densCols(DF$D, DF$MA, colramp = hcl.colors)) +
  theme_bw() +
  xlab(expression(Z-score(Expression~Distance))) +
  ylab(expression(Z-score(Manifold~Distance))) +
  labs(title = 'Expression Distance vs. Manifold Distance', subtitle = parse(text = paste0('rho == ', round(cor(DF$D, DF$MA, method = 'sp'),3)))) +
  theme(plot.title = element_text(face = 2)) + coord_flip()

png('DComparison.png', width = 1500, height = 1500, res = 300)
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
BC <- MASS::boxcox(dE~1)
dE <- dE ^ abs(BC$x[which.max(BC$y)])
zDE <- as.vector(scale(dE))
names(zDE) <- names(dE)

O <- read.csv('results/sym10X500morphineNeuron_Itensor_Dalignment.csv')
Z <- O$Z
names(Z) <- O$gene

tGenes <- intersect(names(zDE), names(Z))

DF <- data.frame(D=(zDE[tGenes]),MA=(Z[tGenes]))

A <- ggplot(DF, mapping = aes(D, MA)) +
  geom_point(col = densCols(DF$D, DF$MA, colramp = hcl.colors)) +
  theme_bw() +
  xlab(expression(Z-score(Tensor~Distance))) +
  ylab(expression(Z-score(Manifold~Distance))) +
  labs(title = 'Expression Distance vs. Manifold Distance', subtitle = parse(text = paste0('rho == ', round(cor(DF$D, DF$MA, method = 'sp'),3)))) +
  theme(plot.title = element_text(face = 2)) + coord_flip()

png('NComparison.png', width = 1500, height = 1500, res = 300)
A
dev.off()
