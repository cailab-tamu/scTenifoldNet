library(Matrix)
library(igraph)
library(ggplot2)
library(patchwork)

O <- list.files('GSD', full.names = TRUE)
O <- paste0(O, '/rankedEdges.csv')

countMatrix <- read.csv('../../inputs/example/GSD/ExpressionData.csv', row.names = 1)
gMean <- rowMeans(countMatrix)


corValue <- sapply(seq_along(O), function(i){
  message(i)
  X <- O[i]
  X <- read.csv(X, sep = '\t')
  colnames(X) <- c('from','to', 'weight')
  X <- graph_from_data_frame(X)
  
  wDegree <- (rowSums(X[]) + colSums(X[]))/2
  #plot(gMean[names(wDegree)], wDegree)
  cor(gMean[names(wDegree)], wDegree, method = 'sp')
})

packageList <- strsplit(O, '/')
packageList <- unlist(lapply(packageList, function(X){X[2]}))

O <- data.frame(P = packageList, corValue = corValue)
O <- O[order(O$corValue),]
O$P <- factor(O$P, levels = O$P)

A <- ggplot(O, aes(corValue, P)) + 
  geom_bar(stat = 'identity') + 
  xlab(expression(rho)) + 
  ylab('Package') + 
  theme_bw() + 
  labs(title = 'GSD') +
  theme(plot.title = element_text(face = 2))

library(scTenifoldNet)
nCells <- seq(100, 2000, 100)
set.seed(1)
cV <- sapply(1:10, function(R){
  sapply(nCells, function(nCell){
    tempCount <- countMatrix[,sample(seq_len(ncol(countMatrix)), size = nCell, replace = TRUE)]
    tCount <- as.matrix(tempCount)
    tempNet <- pcNet(tCount, nComp = 9)
    tDegree <- (colSums(tempNet) + rowSums(tempNet))/2
    tMean <- rowMeans(tCount)      
    tMean <- tMean[names(tDegree)]
    cor(tDegree, tMean, method = 'sp')
  })
})
cV <- t(cV)
cV <- data.frame(Cells = nCells,
           SP = colMeans(cV),
           LB = colMeans(cV) - apply(cV,2,sd),
           UB = colMeans(cV) + apply(cV,2,sd))

B <- ggplot(cV, aes(Cells, SP)) + theme_bw() + 
  geom_point() + 
  geom_errorbar(aes(ymin = LB, ymax=UB)) + 
  ylab(expression(rho)) + 
  xlab('Number of cells') + 
  geom_hline(yintercept = 0, lty = 2, col = 'red') + 
  labs(title = 'SCTENIFOLDNET') +
  theme(plot.title = element_text(face = 2))

png('DegreeExpression.png', width = 800, height = 2000, res = 300)
A / B
dev.off()
