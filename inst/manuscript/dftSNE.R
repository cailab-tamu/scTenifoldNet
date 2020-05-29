library(Seurat)
library(Matrix)
library(ggplot2)
library(scTenifoldNet)
library(RColorBrewer)
library(harmony)

#DF
U <- Read10X('datasets/dermalFibroblasts/unstimulated/')
colnames(U) <- paste0('U_',colnames(U))
U <- scQC(U)
U <- CreateSeuratObject(U)

S <- Read10X('datasets/dermalFibroblasts/stimulated/')
colnames(S) <- paste0('S_', colnames(S))
S <- scQC(S)
S <- CreateSeuratObject(S)

allDF <- merge(U,S)
allDF <- NormalizeData(allDF)
allDF <- ScaleData(allDF)
allDF <- FindVariableFeatures(allDF)
allDF <- RunPCA(allDF, verbose = FALSE)
allDF <- RunHarmony(allDF, group.by.vars = 'orig.ident')
allDF <- RunTSNE(allDF, reduction = 'harmony')
allDF <- RunUMAP(allDF, reduction = 'harmony')
TSNEPlot(allDF)
UMAPPlot(allDF)

gID <- ifelse(Idents(allDF) == 'U', 'Unstimulated', 'Stimulated')
Idents(allDF) <- gID
df <- data.frame(allDF@reductions$tsne@cell.embeddings, Group = as.factor(gID))

png('figures/tsneDF.png', width = 1000, height = 1000, res = 300, pointsize = 30)
ggplot(df, mapping = aes(x = tSNE_1, y = tSNE_2, color = Group)) + geom_point(pch = 20, cex = 1)  + theme_bw() + xlab('t-SNE 1') + ylab('t-SNE 2') + scale_color_brewer(palette='Set1') + theme(legend.position = "none")
dev.off()

VlnPlot(allDF,'RPL3')
VlnPlot(allDF,'RPL10')
VlnPlot(allDF,'RPL32')

png('figures/dfRPL3.png', width = 1000, height = 1000, res = 300)
VlnPlot(allDF, features = 'RPL3', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'RPL3') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/dfRPL10.png', width = 1000, height = 1000, res = 300)
VlnPlot(allDF, features = 'RPL10', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'RPL10') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/dfRPL32.png', width = 1000, height = 1000, res = 300)
VlnPlot(allDF, features = 'RPL32', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'RPL32') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()