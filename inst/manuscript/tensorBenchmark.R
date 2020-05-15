library(scTenifoldNet)
library(ComplexHeatmap)
library(circlize)
library(Matrix)
iMatrix <- readMM('../benchmarking/data/CTL.mtx')
rownames(iMatrix) <- readLines('../benchmarking/data/geneList.txt')

nNet <- 10
nGenes <- nrow(iMatrix)
sGenes <- rownames(iMatrix)
K <- 5
maxError <-  1e-5
maxIter <-  1e3

set.seed(1)
rNetworks <- makeNetworks(iMatrix, nNet = nNet, q = 0.8)

for (i in seq_len(nNet)){
  png(paste0('figures/o',i,'.png'), res = 500, width = 2000, height = 2000, pointsize = 200)
  dataMatrix <- as.matrix(rNetworks[[i]])[1:98,1:98]
  dataMatrix <- (dataMatrix/max(abs(dataMatrix)))
  dataMatrix <- round(dataMatrix,1)
  TP <- sum(dataMatrix[1:40,1:40] > 0) + sum(dataMatrix[41:98,41:98] > 0)
  FP <- sum(dataMatrix[1:40,41:98] > 0) + sum(dataMatrix[41:98, 1:40] > 0)
  TN <- sum(dataMatrix[1:40,41:98] < 0) + sum(dataMatrix[41:98, 1:40] < 0)
  FN <- sum(dataMatrix[1:40,1:40] < 0) + sum(dataMatrix[41:98,41:98] < 0)
  ACC <- round((TP+TN)/(TP+TN+FN+FP),2)
  REC <- round((TP)/((40*40)+(58*58)),2)
  print(ComplexHeatmap::Heatmap(dataMatrix, row_order = 1:98, column_order = 1:98, col= circlize::colorRamp2(breaks = c(-1, 0, 1), 
                                                                                                             colors = c('blue', 'white', 'black')), show_heatmap_legend = FALSE, show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE, border = 'red', 
                                name = 'foo', column_title = paste0('Accuracy = ', ACC, ' | Recall = ', REC)))
  decorate_heatmap_body("foo", {grid.rect(gp = gpar(fill = "transparent", col = "red", lwd = 5))})
  dev.off()
}

tensorX <- array(data = 0, dim = c(nGenes,nGenes,1,nNet))
for(i in seq_len(nNet)){
  tempX <- matrix(0, nGenes, nGenes)
  rownames(tempX) <- colnames(tempX) <- sGenes
  temp <- as.matrix(rNetworks[[i]])
  tGenes <- sGenes[sGenes %in% rownames(temp)]
  tempX[tGenes,tGenes] <- temp[tGenes,tGenes]
  tensorX[,,,i] <- tempX
}

tensorX <- rTensor::as.tensor(tensorX)
set.seed(1)
tensorX <- rTensor::cp(tnsr = tensorX, num_components = K, max_iter = maxIter, tol = maxError)

for (i in seq_len(nNet)){
  png(paste0('figures/t',i,'.png'), res = 500, width = 2000, height = 2000, pointsize = 40)
  dataMatrix <- as.matrix(tensorX$est@data[,,,i])[1:98,1:98]
  dataMatrix <- (dataMatrix/max(abs(dataMatrix)))
  dataMatrix <- round(dataMatrix,1)
  TP <- sum(dataMatrix[1:40,1:40] > 0) + sum(dataMatrix[41:98,41:98] > 0)
  FP <- sum(dataMatrix[1:40,41:98] > 0) + sum(dataMatrix[41:98, 1:40] > 0)
  TN <- sum(dataMatrix[1:40,41:98] < 0) + sum(dataMatrix[41:98, 1:40] < 0)
  FN <- sum(dataMatrix[1:40,1:40] < 0) + sum(dataMatrix[41:98,41:98] < 0)
  ACC <- round((TP+TN)/(TP+TN+FN+FP),2)
  REC <- round((TP)/((40*40)+(58*58)),2)
  print(ComplexHeatmap::Heatmap(dataMatrix, row_order = 1:98, column_order = 1:98, col= circlize::colorRamp2(breaks = c(-1, 0, 1), 
        colors = c('blue', 'white', 'black')), show_heatmap_legend = FALSE, show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE, border = 'forestgreen', 
        name = 'foo', column_title = paste0('Accuracy = ', ACC, ' | Recall = ', REC)))
  decorate_heatmap_body("foo", {grid.rect(gp = gpar(fill = "transparent", col = "forestgreen", lwd = 5))})
  dev.off()
}

outO <- (rNetworks[[1]] + rNetworks[[2]] + rNetworks[[3]] + rNetworks[[4]] + rNetworks[[5]] + rNetworks[[6]] + rNetworks[[7]] + rNetworks[[8]] + rNetworks[[9]] + rNetworks[[10]])/10
outO <- outO/max(abs(outO))
outO <- round(outO,1)
png('figures/avgO.png', res = 500, width = 2000, height = 2000, pointsize = 40)
dataMatrix <- as.matrix(outO)[1:98,1:98]
TP <- sum(dataMatrix[1:40,1:40] > 0) + sum(dataMatrix[41:98,41:98] > 0)
FP <- sum(dataMatrix[1:40,41:98] > 0) + sum(dataMatrix[41:98, 1:40] > 0)
TN <- sum(dataMatrix[1:40,41:98] < 0) + sum(dataMatrix[41:98, 1:40] < 0)
FN <- sum(dataMatrix[1:40,1:40] < 0) + sum(dataMatrix[41:98,41:98] < 0)
ACC <- round((TP+TN)/(TP+TN+FN+FP),2)
REC <- round((TP)/((40*40)+(58*58)),2)
print(ComplexHeatmap::Heatmap(dataMatrix, row_order = 1:98, column_order = 1:98, col= circlize::colorRamp2(breaks = c(-1, 0, 1), 
                                                                                                           colors = c('blue', 'white', 'black')), show_heatmap_legend = FALSE, show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE, border = 'red', 
                              name = 'foo', column_title = paste0('Accuracy = ', ACC, ' | Recall = ', REC)))
decorate_heatmap_body("foo", {grid.rect(gp = gpar(fill = "transparent", col = "red", lwd = 5))})
dev.off()

outT <- (tensorX$est@data[,,,1] + tensorX$est@data[,,,2] + tensorX$est@data[,,,3] + tensorX$est@data[,,,4] + tensorX$est@data[,,,5] + tensorX$est@data[,,,6] + tensorX$est@data[,,,7] + tensorX$est@data[,,,8] + tensorX$est@data[,,,9] + tensorX$est@data[,,,10])/10

outT <- outT/max(abs(outT))
outT <- round(outT,1)

png('figures/avgT.png', res = 500, width = 2000, height = 2000, pointsize = 40)
dataMatrix <- as.matrix(outT)[1:98,1:98]
TP <- sum(dataMatrix[1:40,1:40] > 0) + sum(dataMatrix[41:98,41:98] > 0)
FP <- sum(dataMatrix[1:40,41:98] > 0) + sum(dataMatrix[41:98, 1:40] > 0)
TN <- sum(dataMatrix[1:40,41:98] < 0) + sum(dataMatrix[41:98, 1:40] < 0)
FN <- sum(dataMatrix[1:40,1:40] < 0) + sum(dataMatrix[41:98,41:98] < 0)
ACC <- round((TP+TN)/(TP+TN+FN+FP),2)
REC <- round((TP)/((40*40)+(58*58)),2)
print(ComplexHeatmap::Heatmap(dataMatrix, row_order = 1:98, column_order = 1:98, col= circlize::colorRamp2(breaks = c(-1, 0, 1), 
                                                                                                           colors = c('blue', 'white', 'black')), show_heatmap_legend = FALSE, show_row_names = FALSE, show_column_names = FALSE, use_raster = TRUE, border = 'forestgreen', 
                              name = 'foo', column_title = paste0('Accuracy = ', ACC, ' | Recall = ', REC)))
decorate_heatmap_body("foo", {grid.rect(gp = gpar(fill = "transparent", col = "forestgreen", lwd = 5))})
dev.off()
