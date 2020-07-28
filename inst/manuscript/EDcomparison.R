library('Matrix')
library(igraph)
library(ggplot2)
library(patchwork)

eX <- readMM('datasets/morphineNeurons/mock/control.mtx')
rownames(eX) <- readLines('datasets/morphineNeurons/mock/controlGenes.tsv')
eX <- rowMeans(eX)

eY <- readMM('datasets/morphineNeurons/morphine/morphine.mtx')
rownames(eY) <- readLines('datasets/morphineNeurons/morphine/morphineGenes.tsv')
eY <- rowMeans(eY)

nX <- readMM('results/tensorOutput/X_10X500morphineNeuron_Itensor.mtx')
rownames(nX) <- readLines('results/tensorOutput/genes_10X500morphineNeuron_Itensor.mtx')
colnames(nX) <- readLines('results/tensorOutput/genes_10X500morphineNeuron_Itensor.mtx')
nX <- graph_from_adjacency_matrix(nX!=0)
nX <- degree(nX)

nY <- readMM('results/tensorOutput/Y_10X500morphineNeuron_Itensor.mtx')
rownames(nY) <- readLines('results/tensorOutput/genes_10X500morphineNeuron_Itensor.mtx')
colnames(nY) <- readLines('results/tensorOutput/genes_10X500morphineNeuron_Itensor.mtx')
nY <- graph_from_adjacency_matrix(nY!=0)
nY <- degree(nY)

O <- read.csv('results/sym10X500morphineNeuron_Itensor_Dalignment.csv')
Z <- O$distance
names(Z) <- O$gene
sharedGenes <- intersect(intersect(intersect(names(eX), names(nX)), intersect(names(eY), names(nY))), names(Z))


DF <- data.frame(eX = eX[sharedGenes], eY = eY[sharedGenes], nX = nX[sharedGenes], nY = nY[sharedGenes], Z = Z[sharedGenes])
DF <- DF[(DF[,3] != 0) & (DF[,4] != 0),]
DF$eX <- log1p(DF$eX)
DF$eY <- log1p(DF$eY)
DF$nX <- log10(DF$nX)
DF$nY <- log10(DF$nY)
DF$Z <- log10(DF$Z)

A <- ggplot(DF, mapping = aes(eX, nX)) + 
  geom_point(col = densCols(DF$eX, DF$nX, colramp = hcl.colors), cex = 0.5) + 
  theme_bw() + 
  xlab(expression(log(Average~Expression + 1))) + 
  ylab(expression(log[10](Degree))) + 
  labs(title = 'Control', subtitle = parse(text = paste0('rho == ', round(cor(DF$eX, DF$nX, method = 'sp'),3)))) + 
  theme(plot.title = element_text(face = 2))
B <- ggplot(DF, mapping = aes(eY, nY)) + 
  geom_point(col = densCols(DF$eY, DF$nY, colramp = hcl.colors), cex = 0.5) + 
  theme_bw() + 
  xlab(expression(log(Average~Expression + 1))) + 
  ylab(expression(log[10](Degree))) + 
  labs(title = 'Morphine', subtitle = parse(text = paste0('rho == ', round(cor(DF$eY, DF$nY, method = 'sp'),3)))) + 
  theme(plot.title = element_text(face = 2))
C <- ggplot(DF, mapping = aes(eX, Z)) + 
  geom_point(col = densCols(DF$eX, DF$Z, colramp = hcl.colors), cex = 0.5) + 
  theme_bw() + 
  labs(subtitle = parse(text = paste0('rho == ', round(cor(DF$eX, DF$Z, method = 'sp'),3)))) +
  xlab(expression(log(Average~Expression + 1))) + 
  ylab(expression(log[10](Distance)))
D <- ggplot(DF, mapping = aes(eY, Z)) + geom_point(col = densCols(DF$eY, DF$Z, colramp = hcl.colors), cex = 0.5) + 
  theme_bw() + 
  xlab(expression(log(Average~Expression + 1))) + 
  ylab(expression(log[10](Distance)))

png('EDComparison.png', width = 1500, height = 1500, res = 300)
A + B + C + D
dev.off()
