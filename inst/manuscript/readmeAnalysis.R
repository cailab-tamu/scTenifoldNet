library(scTenifoldNet)

nCells = 2000
nGenes = 100
set.seed(1)
X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
X <- round(X)
X <- matrix(X, ncol = nCells)
mean(X != 0)
rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))

Y <- X
Y[10,] <- Y[50,]
Y[2,] <- Y[11,]
Y[3,] <- Y[5,]

outputH0 <- scTenifoldNet(X = X, Y = X,
                          nc_nNet = 10, nc_nCells = 500,
                          td_K = 3, qc_minLibSize = 30,
                          dc_minFC = 0)

outputHA <- scTenifoldNet(X = X, Y = Y,
                          nc_nNet = 10, nc_nCells = 500,
                          td_K = 3, qc_minLibSize = 30,
                          dc_minFC = 0)

png('figures/readmeExample.png', width = 900, height = 900, res = 300, pointsize = 10)
par(mfrow = c(1,2), mar=c(2.5,2.5,1,1), mgp = c(1.5,0.5,0))
plotQQ <- function(D, ...){
  P <- D$diffRegulation$p.adj
  gColor <- ifelse(P < 0.05, 'red', 'black')
  qqnorm(D$diffRegulation$Z, pch = 20, col= gColor, ...)
  qqline(D$diffRegulation$Z)
  legend('bottomright', legend = 'FDR < 0.05', col = 'red', pch = 20, bty = 'n', cex = 0.7)
}
plotQQ(outputH0, main = 'Original vs Original', xlab = 'Normal Quantiles')
plotQQ(outputHA, main = 'Original vs Perturbed', xlab = 'Normal Quantiles')
dev.off()
