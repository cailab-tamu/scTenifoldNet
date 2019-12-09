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

oSample <- CreateSeuratObject(oSample)
ySample <- CreateSeuratObject(ySample)

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

VlnPlot(allAging,'Meis2')
VlnPlot(allAging,'Celf2')
VlnPlot(allAging,'mt-Cytb')

png('figures/agingMEIS2.png', width = 1000, height = 1000, res = 300)
VlnPlot(allAging, features = 'Meis2', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'MEIS2') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/agingCELF2.png', width = 1000, height = 1000, res = 300)
VlnPlot(allAging, features = 'Celf2', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'CELF2') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/agingCOX7C.png', width = 1000, height = 1000, res = 300)
VlnPlot(allAging, features = 'Cox7c', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'COX7C') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()
