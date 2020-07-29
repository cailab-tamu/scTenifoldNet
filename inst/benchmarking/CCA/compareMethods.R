library(Matrix)
library(GGally)

load('CCA_morphineNeuron.Rdata')
MA <- read.csv('sym10X500morphineNeuron_Itensor_Dalignment.csv', row.names = 1)

zMA <- MA$Z
names(zMA) <- MA$gene

zCCA2 <- dRCCA_d2$Z
names(zCCA2) <- dRCCA_d2$gene

zCCA3 <- dRCCA_d3$Z
names(zCCA3) <- dRCCA_d3$gene

zCCA4 <- dRCCA_d4$Z
names(zCCA4) <- dRCCA_d4$gene

zCCA5 <- dRCCA_d5$Z
names(zCCA5) <- dRCCA_d5$gene

sGenes <- intersect(names(zMA), names(zCCA2))
sGenes <- intersect(sGenes, names(zCCA3))
sGenes <- intersect(sGenes, names(zCCA4))

ALL <- data.frame(MA=zMA[sGenes], CCA2 = zCCA2[sGenes], CCA3 = zCCA3[sGenes], CCA4 = zCCA4[sGenes], CCA5 = zCCA5[sGenes])
rALL <- as.data.frame(apply(ALL,2,rank))

png('MA_CCA.png', width = 2000, height = 2000, res = 300)
ggpairs(rALL)
dev.off()
