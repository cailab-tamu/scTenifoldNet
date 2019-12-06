library(Seurat)
library(Matrix)
library(ggplot2)
library(scTenifoldNet)
library(RColorBrewer)

# Aging
ySample <- readMM('datasets/agingNeurons/Y/yMatrix.mtx')
rownames(ySample) <- readLines('datasets/agingNeurons/Y/yGenes.txt')
colnames(ySample) <- paste0('Y_',readLines('datasets/agingNeurons/Y/yBarcodes.txt'))
ySample <- scQC(ySample)

oSample <- readMM('datasets/agingNeurons/O/oMatrix.mtx')
rownames(oSample) <- readLines('datasets/agingNeurons/O/oGenes.txt')
colnames(oSample) <- paste0('O_',readLines('datasets/agingNeurons/O/oBarcodes.txt'))
oSample <- scQC(oSample)

#par(mfrow=c(1,2))
par(mar=c(3,3,1,1), mgp = c(1.5,0.5,0))
dta <- cbind(ySample['Malat1',], ySample['Ggt7',])
dta[dta == 0] <- NA
plot(dta, col = 'black', pch = 20, cex = 0.7, xlab = 'MALAT1', ylab = 'GGT7')
abline(lm(dta[,2]~dta[,1]))
cor(dta[,1],dta[,2], use = 'com', method = 'sp')

dta <- cbind(oSample['Malat1',], oSample['Ggt7',])
dta[dta == 0] <- NA
points(dta, col = 'red', pch = 20, cex = 0.7)
abline(lm(dta[,2]~dta[,1]), col = 'red')
cor(dta[,1],dta[,2], use = 'com', method = 'sp')

dta <- cbind(c(ySample['Malat1',],oSample['Malat1',]), c(ySample['Ggt7',],oSample['Ggt7',]), c(rep(0, ncol(ySample)), rep(1, ncol(oSample))))
colnames(dta) <- c('X', 'Y', 'G')
dta <- as.data.frame(dta)
dta <- dta[rowSums(dta[,1:2] == 0) == 0,]
dta$G <- as.factor(dta$G)
plot(dta[,1:2])
stargazer(lme4::lmer(Y~X+(1|G), dta), type = 'text')

plot(dta[,1:2])

par(mar=c(3,3,1,1), mgp = c(1.5,0.5,0))
dta <- cbind(ySample['Gad1',], ySample['Stau2',])
dta[dta == 0] <- NA
plot(dta, col = 'black', pch = 20, cex = 0.7, xlab = 'MALAT1', ylab = 'GGT7')
abline(lm(dta[,2]~dta[,1]))
cor(dta[,1],dta[,2], use = 'com', method = 'sp')

dta <- cbind(oSample['Gad1',], oSample['Stau2',])
dta[dta == 0] <- NA
points(dta, col = 'red', pch = 20, cex = 0.7)
abline(lm(dta[,2]~dta[,1]), col = 'red')
cor(dta[,1],dta[,2], use = 'com', method = 'sp')




oSample <- CreateSeuratObject(oSample)

allAging <- merge(ySample, oSample)
allAging <- NormalizeData(allAging)
allAging <- ScaleData(allAging)
allAging <- FindVariableFeatures(allAging)
allAging <- RunPCA(allAging, verbose = FALSE)
allAging <- RunTSNE(allAging)


gID <- ifelse(Idents(allAging) == 'O', 'Old', 'Young')
Idents(allAging) <- gID
df <- data.frame(allAging@reductions$tsne@cell.embeddings, Group = as.factor(gID))

png('figures/tsneAging.png', width = 1000, height = 1000, res = 300, pointsize = 30)
ggplot(df, mapping = aes(x = tSNE_1, y = tSNE_2, color = Group)) + geom_point(pch = 20, cex = 1)  + theme_bw() + xlab('t-SNE 1') + ylab('t-SNE 2') + scale_color_brewer(palette='Set1') + theme(legend.position = "none")
dev.off()

VlnPlot(allAging,'mt-Co2')
VlnPlot(allAging,'Celf2')
VlnPlot(allAging,'mt-Cytb')

png('figures/agingCO2.png', width = 1000, height = 1000, res = 300)
VlnPlot(allAging, features = 'mt-Co2', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'MT-CO2') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/agingCELF2.png', width = 1000, height = 1000, res = 300)
VlnPlot(allAging, features = 'Celf2', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'CELF2') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/agingCYTB.png', width = 1000, height = 1000, res = 300)
VlnPlot(allAging, features = 'mt-Cytb', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'MT-CYTB') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()
