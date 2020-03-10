library(Seurat)
library(Matrix)
library(ggplot2)
library(scTenifoldNet)
library(RColorBrewer)

# Morphine
morphineC <- readMM('../manuscript/datasets/morphineNeurons/mock/control.mtx')
colnames(morphineC) <- readLines('datasets/morphineNeurons/mock/controlBarcodes.tsv')
rownames(morphineC) <- readLines('datasets/morphineNeurons/mock/controlGenes.tsv')
morphineC <- scQC(morphineC)
morphineC <- CreateSeuratObject(morphineC)

morphineM <- readMM('../manuscript/datasets/morphineNeurons/morphine/morphine.mtx')
colnames(morphineM) <- readLines('datasets/morphineNeurons/morphine/morphineBarcodes.tsv')
rownames(morphineM) <- readLines('datasets/morphineNeurons/morphine/morphineGenes.tsv')
morphineM <- scQC(morphineM)
morphineM <- CreateSeuratObject(morphineM)

allMorphine <- merge(morphineC, morphineM)
allMorphine <- NormalizeData(allMorphine)
allMorphine <- ScaleData(allMorphine)
allMorphine <- FindVariableFeatures(allMorphine)
allMorphine <- RunPCA(allMorphine, verbose = FALSE)
allMorphine <- RunTSNE(allMorphine)
allMorphine <- RunUMAP(allMorphine, dims = 1:50)
UMAPPlot(allMorphine)

gID <- gsub('[[:digit:]]','',Idents(allMorphine))
gID <- ifelse(gID == 'C', 'Control', 'Morphine')
df <- data.frame(allMorphine@reductions$tsne@cell.embeddings, Group = as.factor(gID))

png('figures/tsneMorphine.png', width = 1000, height = 1000, res = 300, pointsize = 30)
ggplot(df, mapping = aes(x = tSNE_1, y = tSNE_2, color = Group)) + geom_point(pch = 20, cex = 1)  + theme_bw() + xlab('t-SNE 1') + ylab('t-SNE 2') + scale_color_brewer(palette='Set1') + theme(legend.position = "none")
dev.off()

png('figures/morphineADCY5.png', width = 1000, height = 1000, res = 300)
VlnPlot(allMorphine, features = 'Adcy5', cols = c(rep('#E41A1C',4), rep("#377EB8",4)), pt.size = 0.25)+ theme_bw() + labs(title = 'Adcy5') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/morphineGABRG1.png', width = 1000, height = 1000, res = 300)
VlnPlot(allMorphine, features = 'Gabrg1', cols = c(rep('#E41A1C',4), rep("#377EB8",4)), pt.size = 0.25)+ theme_bw() + labs(title = 'Gabrg1') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/morphinePDE1B.png', width = 1000, height = 1000, res = 300)
VlnPlot(allMorphine, features = 'Pde1b', cols = c(rep('#E41A1C',4), rep("#377EB8",4)), pt.size = 0.25)+ theme_bw() + labs(title = 'Pde1b') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()
