library(Matrix)
library(scTenifoldNet)
library(GGally)
library(UpSetR)
library(patchwork)
library(ggrepel)

load('DM100.RData')
J100 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n100 <- nrow(O$diffRegulation)
z100 <- O$diffRegulation$Z
names(z100) <- O$diffRegulation$gene

load('DM200.RData')
J200 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n200 <- nrow(O$diffRegulation)
z200 <- O$diffRegulation$Z
names(z200) <- O$diffRegulation$gene

load('DM2000.RData')
J2000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n2000 <- nrow(O$diffRegulation)
z2000 <- O$diffRegulation$Z
names(z2000) <- O$diffRegulation$gene

load('DM3000.RData')
J3000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n3000 <- nrow(O$diffRegulation)
z3000 <- O$diffRegulation$Z
names(z3000) <- O$diffRegulation$gene

load('DM4000.RData')
J4000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n4000 <- nrow(O$diffRegulation)
z4000 <- O$diffRegulation$Z
names(z4000) <- O$diffRegulation$gene

load('DM5000.RData')
J5000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n5000 <- nrow(O$diffRegulation)
z5000 <- O$diffRegulation$Z
names(z5000) <- O$diffRegulation$gene

load('DM6000.RData')
J6000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n6000 <- nrow(O$diffRegulation)
z6000 <- O$diffRegulation$Z
names(z6000) <- O$diffRegulation$gene

load('DM7000.RData')
J7000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000 <- nrow(O$diffRegulation)
z7000 <- O$diffRegulation$Z
names(z7000) <- O$diffRegulation$gene

sGenes <- intersect(names(z100), names(z200))
sGenes <- intersect(sGenes, names(z2000))
sGenes <- intersect(sGenes, names(z3000))
sGenes <- intersect(sGenes, names(z4000))
sGenes <- intersect(sGenes, names(z5000))
sGenes <- intersect(sGenes, names(z6000))
sGenes <- intersect(sGenes, names(z7000))

allZ <- data.frame(c100 = z100[sGenes], c200 = z200[sGenes], c2000 = z2000[sGenes], c3000 = z3000[sGenes], c4000 = z4000[sGenes], c5000 = z5000[sGenes], c6000 = z6000[sGenes], c7000 = z7000[sGenes])

PCA <- prcomp(t(allZ))$x
PCA <- as.data.frame(PCA)
PCA$Cells <- gsub('c','',rownames(PCA))
A <- ggplot(PCA, aes(PC1,PC2, label=Cells)) + geom_point() + geom_text_repel() + theme_bw()


cValues <- cor(allZ, method = 'sp')
diag(cValues) <- NA
rownames(cValues) <- colnames(cValues) <- c('100', '200', '2000', '3000', '4000', '5000', '6000', '7000')
cValues <- reshape2::melt(cValues)
cValues$Var1 <- as.factor(cValues$Var1)
cValues$Var2 <- as.factor(cValues$Var2)
colnames(cValues)[3] <- 'SCC'
B <- ggplot(cValues, aes(Var1, Var2, fill=SCC)) + geom_tile() + theme_minimal() + xlab('Number of Cells') + ylab('Number of Cells')


allJ <- fromList(list(c100 = J100, c200 = J200, c2000 = J2000, c3000 = J3000, c4000 = J4000, c5000 = J5000, c6000 = J6000, c7000 = J7000))
jValues <- sapply(seq_len(ncol(allJ)), function(i){
  sapply(seq_len(ncol(allJ)), function(j){
    tempX <- cbind(allJ[,i], allJ[,j])
    tempX <- rowSums(tempX)
    sum(tempX == 2)/sum(tempX > 0)
  })
})
#jValues[lower.tri(jValues, diag = TRUE)] <- NA
diag(jValues) <- NA
rownames(jValues) <- colnames(jValues) <- c('100', '200', '2000', '3000', '4000', '5000', '6000', '7000')
jValues <- reshape2::melt(jValues)
jValues$Var1 <- as.factor(jValues$Var1)
jValues$Var2 <- as.factor(jValues$Var2)
colnames(jValues)[3] <- 'Jaccard'
C <- ggplot(jValues, aes(Var1, Var2, fill=Jaccard)) + geom_tile() + theme_minimal() + xlab('Number of Cells') + ylab('Number of Cells')


nGenes <- rbind(n100, n200, n2000, n3000, n4000, n5000, n6000, n7000)
nGenes <- as.data.frame(nGenes)
nGenes$Cells <- c('100', '200', '2000', '3000', '4000', '5000', '6000', '7000')
D <- ggplot(nGenes, aes(y = Cells, x = V1)) + geom_bar(stat = 'identity') + theme_bw() + xlab('Number of Genes') + ylab('Number of Cells')  + labs(title = 'Mouse Neurons') + theme(plot.title = element_text(face = 2))

png('DownSampling.png', width = 2400, height = 2000, res = 300)
D + B + A + C
dev.off()
