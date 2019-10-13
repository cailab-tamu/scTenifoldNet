library(Matrix)

Accuracy <- function(X){
  c1p <- sum(X[1:40,1:40] > 0)
  c2p <- sum(X[41:98,41:98] > 0)
  c1n <- sum(X[1:40,41:98] < 0)
  c2n <- sum(X[41:98,1:40] < 0)
  (c1p+c1n+c2p+c2n)/sum(X[1:98,1:98] != 0)
}

reCall <- function(X){
  c1p <- sum(X[1:40,1:40] > 0)
  c2p <- sum(X[41:98,41:98] > 0)
  (c1p+c2p)/((40*40)+(58*58))
}

nCells <- c(500, 1000, 2000, 3000)
for(c in seq_along(nCells)){
  png(paste0('figures/tsrPR_',nCells[c],'cells.png'), width = 900*4, height = 900*2, res = 300)
  layout(matrix(c(1,2,3,4,1,2,3,4,1,2,3,4,5,6,7,8,5,6,7,8,5,6,7,8,9,9,9,9),byrow = TRUE, ncol = 4))
  par(mar=c(3,3,1,1), mgp=c(1.5,0.5,0))
  netColors <- RColorBrewer::brewer.pal(6, 'Paired')
  nThresholds <- c(0.95,0.9, 0.85,0.8)
  nNetworks <- c(5,10,15,20,50,100)
  for(j in seq_along(nThresholds)){
    plot(0, type = 'n', xlim = c(0,1), ylim = c(0.6,1), xlab = 'Recall', ylab = 'Precision', main = nThresholds[j])
    for(i in seq_along(nNetworks)){
      X <- readMM(paste0('networks/TENSOR/tsr_',nCells[c],'cells_',nNetworks[i],'nets_',nThresholds[j],'q.mtx'))
      X <- as.matrix(X)
      out <- t(sapply(seq(1,0,-0.05), function(Z){
        X[abs(X) < quantile(abs(X),Z)] <- 0
        c(reCall(X),Accuracy(X))
      }))
      points(out, type = 'b', col = netColors[i], pch = (15+i))
    }
  }
  for(j in seq_along(nThresholds)){
    plot(0, type = 'n', xlim = c(0,1), ylim = c(0.6,1), xlab = 'Recall', ylab = 'Precision', main = nThresholds[j])
    for(i in seq_along(nNetworks)){
      X <- readMM(paste0('networks/TENSOR/avg_',nCells[c],'cells_',nNetworks[i],'nets_',nThresholds[j],'q.mtx'))
      X <- as.matrix(X)
      out <- t(sapply(seq(1,0,-0.05), function(Z){
        X[abs(X) < quantile(abs(X),Z)] <- 0
        c(reCall(X),Accuracy(X))
      }))
      points(out, type = 'b', col = netColors[i], pch = (15+i))
    }
  }
  plot.new()
  par(mar=c(0,0,0,0))
  legend('top', legend = paste0(nNetworks, ' Networks'), pch=c(16:22), bty = 'n', ncol = 6, col = netColors)
  # plot.new()
  # par(mar=c(0,0,2,0))
  # legend('top',legend = '', bty = 'n', cex = 2, title = paste0(nCells[c], ' CELLS'))
  dev.off()
}


library(Matrix)
library(ComplexHeatmap)
AVG <- readMM('networks/TENSOR/avg_1000cells_10nets_0.85q.mtx')
TSR <- readMM('networks/TENSOR/tsr_1000cells_10nets_0.85q.mtx')

AVG <- as.matrix(AVG)
AVG <- AVG/max(abs(AVG))
TSR <- as.matrix(TSR)
AVG[AVG < 0.1] <- -1
TSR[TSR < 0.1] <- -1
library(circlize)
Heatmap(AVG, column_order=1:100, row_order = 1:100, show_heatmap_legend = FALSE, col = colorRamp2(c(-0.1, 0, 1), c("blue", "white", "red")))+
  Heatmap(TSR, column_order=1:100, row_order = 1:100, name = 'PCr', col=colorRamp2(c(-1, 0, 1), c("blue", "white", "red")))

