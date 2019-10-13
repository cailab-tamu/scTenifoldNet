library(Matrix)
library(microbenchmark)
library(scTenifoldNet)

countMatrix <- readMM('data/CTL.mtx')
countMatrix <- as.matrix(countMatrix)
rownames(countMatrix) <- readLines('data/geneList.txt')


PCr <- function(countMatrix){
  countMatrix <- countMatrix[,sample(seq_len(ncol(countMatrix)),1000)]
  scTenifoldNet::pcNet(countMatrix)
}

SCC <- function(countMatrix){
  countMatrix <- countMatrix[,sample(seq_len(ncol(countMatrix)),1000)]
  stats::cor(countMatrix, method = 'spearman')
}

MI <- function(countMatrix){
  countMatrix <- countMatrix[,sample(seq_len(ncol(countMatrix)),1000)]
  minet::minet(t(countMatrix))
}

GENIE3 <- function(countMatrix){
  countMatrix <- countMatrix[,sample(seq_len(ncol(countMatrix)),1000)]
  GENIE3::GENIE3(countMatrix)
}

timeBenchmark <- microbenchmark(
PCr(countMatrix),
SCC(countMatrix),
MI(countMatrix),
GENIE3(countMatrix),times = 10, unit = 's')

timeBenchmark <- summary(timeBenchmark)
write.csv(timeBenchmark, row.names = FALSE, file = 'results/timeBenchmark.csv')
barColor <- RColorBrewer::brewer.pal(4,'Set1')
orderMethod <- order(timeBenchmark$mean)
nameMethod <- c('PCr', 'SCC', 'MI', 'GENIE3')
barPlot <- timeBenchmark$mean[orderMethod]
names(barPlot) <- nameMethod[orderMethod]

png('figures/timeBenchmarking.png', width = 2400, height = 500, res = 300)
par(mar=c(3,5,1,1), mgp=c(1.5,0.5,0))
barplot(barPlot, horiz = TRUE, xlim = c(0, 1.2*max(timeBenchmark$max)), xlab = 'Seconds', las = 1, col= barColor)
dev.off()
