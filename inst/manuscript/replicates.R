library(enrichR)
library(scTenifoldNet)

par(mfrow=c(1,3))
A <- read.csv('results/10X1000DF_Itensor_Dalignment.csv', row.names = 1)[,1:30]
B <- read.csv('results/10X500DF_Itensor_Dalignment.csv', row.names = 1)[,1:30]
A <- dCoexpression(A, minDist = 1e-5)
B <- dCoexpression(B, minDist = 1e-5)
rownames(A) <- A$gene
rownames(B) <- B$gene
sharedGenes <- intersect(A$gene, B$gene)
A <- A[sharedGenes,]
B <- B[sharedGenes,]
genePCH <- ifelse(A$p.value < 0.05 & B$p.value < 0.05, 8, 20)
geneColors <- densCols(cbind(A$Z,B$Z))
geneColors[A$p.value < 0.05 | B$p.value < 0.05] <- 'red'
par(mar=c(3,3,1,1), mgp=c(1.5,0.5,0))
plot(A$Z, B$Z, col= geneColors, pch = genePCH, xlab = 'Z-score(Distance) 1000 Cells', ylab = 'Z-score(Distance) 500 Cells', xlim=c(0,15), ylim = c(0,15))
cor.test(log10(1+A$distance), log10(1+B$distance), method = 'sp')
abline(0,1, lty = 2)
abline(v = min(A$Z[A$p.value < 0.05]), lty = 3, col = 'red')
abline(h = min(B$Z[B$p.value < 0.05]), lty = 3, col = 'red')


A <- read.csv('results/10X1000ANEURON_Itensor_Dalignment.csv', row.names = 1)[,1:30]
B <- read.csv('results/10X500ANEURON_Itensor_Dalignment.csv', row.names = 1)[,1:30]
A <- dCoexpression(A, minDist = 1e-5)
B <- dCoexpression(B, minDist = 1e-5)
rownames(A) <- A$gene
rownames(B) <- B$gene
sharedGenes <- intersect(A$gene, B$gene)
A <- A[sharedGenes,]
B <- B[sharedGenes,]
genePCH <- ifelse(A$p.value < 0.05 & B$p.value < 0.05, 8, 20)
geneColors <- densCols(cbind(A$Z,B$Z))
geneColors[A$p.value < 0.05 | B$p.value < 0.05] <- 'red'
par(mar=c(3,3,1,1), mgp=c(1.5,0.5,0))
plot(A$Z, B$Z, col= geneColors, pch = genePCH, xlab = 'Z-score(Distance) 1000 Cells', ylab = 'Z-score(Distance) 500 Cells', xlim=c(0,15), ylim = c(0,15))
cor.test(log10(1+A$distance), log10(1+B$distance), method = 'sp')
abline(0,1, lty = 2)
abline(v = min(A$Z[A$p.value < 0.05]), lty = 3, col = 'red')
abline(h = min(B$Z[B$p.value < 0.05]), lty = 3, col = 'red')


A <- read.csv('results/10X1000morphineNeuron_Itensor_Dalignment.csv', row.names = 1)[,1:30]
B <- read.csv('results/10X500morphineNeuron_Itensor_Dalignment.csv', row.names = 1)[,1:30]
A <- dCoexpression(A, minDist = 1e-5)
B <- dCoexpression(B, minDist = 1e-5)
rownames(A) <- A$gene
rownames(B) <- B$gene
sharedGenes <- intersect(A$gene, B$gene)
A <- A[sharedGenes,]
B <- B[sharedGenes,]
genePCH <- ifelse(A$p.value < 0.05 & B$p.value < 0.05, 8, 20)
geneColors <- densCols(cbind(A$Z,B$Z))
geneColors[A$p.value < 0.05 | B$p.value < 0.05] <- 'red'
par(mar=c(3,3,1,1), mgp=c(1.5,0.5,0))
plot(A$Z, B$Z, col= geneColors, pch = genePCH, xlab = 'Z-score(Distance) 1000 Cells', ylab = 'Z-score(Distance) 500 Cells', xlim=c(0,15), ylim = c(0,15))
cor.test(log10(1+A$distance), log10(1+B$distance), method = 'sp')
abline(0,1, lty = 2)
abline(v = min(A$Z[A$p.value < 0.05]), lty = 3, col = 'red')
abline(h = min(B$Z[B$p.value < 0.05]), lty = 3, col = 'red')
