library(devtools)
load_all('.')

nCells = 2000
nGenes = 100
set.seed(1)
X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
X <- round(X)
X <- matrix(X, ncol = nCells)
rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))

Y <- X
Y[10,] <- Y[50,]
Y[2,] <- Y[11,]
Y[3,] <- Y[5,]

outputHA_td3 <- scTenifold3Net(X = X, Y = Y,
                          nc_nNet = 10, nc_nCells = 500,
                          td_K = 3, qc_minLibSize = 30)

outputHA_td4 <- scTenifoldNet(X = X, Y = Y,
                            nc_nNet = 10, nc_nCells = 500,
                            td_K = 3, qc_minLibSize = 30)

D1 <- outputHA_td3$diffRegulation$distance
names(D1) <- outputHA_td3$diffRegulation$gene

D2 <- outputHA_td4$diffRegulation$distance
names(D2) <- outputHA_td4$diffRegulation$gene

D2 <- D2[names(D1)]

png('tensorDimension.png', width = 1500, height = 1500, res = 300)
par(mar=c(3,3,1,1), mgp=c(1.5,0.5,0))
gCol <- ifelse(outputHA_td3$diffRegulation$p.adj < 0.1 | outputHA_td4$diffRegulation$p.adj < 0.1, 'red', 'black')
plot(D2,D1, xlab = 'Tensor 4D', ylab = 'Tensor 3D', pch = 16, col = gCol)
abline(0,1, col = 'red', lty = 2)
dev.off()
