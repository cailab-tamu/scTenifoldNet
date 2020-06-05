library(Seurat)
library(Matrix)
library(ggplot2)
library(scTenifoldNet)
library(RColorBrewer)
library(harmony)

#NKX21
WT <- readMM('datasets/NKX2-1/wtAT1.mtx')
rownames(WT) <- readLines('datasets/NKX2-1/genes_wtAT1.txt')
colnames(WT) <- readLines('datasets/NKX2-1/barcodes_wtAT1.txt')
WT <- CreateSeuratObject(WT, project = 'WT')

KO <- readMM('datasets/NKX2-1/koAT1.mtx')
colnames(KO) <- readLines('datasets/NKX2-1/barcodes_koAT1.txt')
rownames(KO) <- readLines('datasets/NKX2-1/genes_koAT1.txt')
KO <- CreateSeuratObject(KO, project = 'KO')

ALL <- merge(WT,KO)
ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)
ALL <- FindVariableFeatures(ALL)
ALL <- RunPCA(ALL, verbose = FALSE)
ALL <- RunHarmony(ALL, group.by.vars = 'orig.ident')
ALL <- RunTSNE(ALL, reduction = 'harmony')
ALL <- RunUMAP(ALL, reduction = 'harmony', dims = 1:20)
TSNEPlot(ALL)
UMAPPlot(ALL)

gID <- Idents(ALL)
df <- data.frame(ALL@reductions$tsne@cell.embeddings, Group = as.factor(gID))

png('figures/tsneNKX21.png', width = 1000, height = 1000, res = 300, pointsize = 30)
ggplot(df, mapping = aes(x = tSNE_1, y = tSNE_2, color = Group)) + geom_point(pch = 20, cex = 1)  + theme_bw() + xlab('t-SNE 1') + ylab('t-SNE 2') + scale_color_brewer(palette='Set1') + theme(legend.position = "none")
dev.off()

VlnPlot(ALL,'Fxyd3')
VlnPlot(ALL,'Cd24a')
VlnPlot(ALL,'Fau')

png('figures/nkx21FAU.png', width = 1000, height = 1000, res = 300)
VlnPlot(ALL, features = 'Fau', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'Fau') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/nkx21CD24A.png', width = 1000, height = 1000, res = 300)
VlnPlot(ALL, features = 'Cd24a', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'Cd24a') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/nkx21EEF1A1.png', width = 1000, height = 1000, res = 300)
VlnPlot(ALL, features = 'Eef1a1', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'Eef1a1') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()