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

VlnPlot(allDF,'SOD2')
VlnPlot(allDF,'NFKBIA')
VlnPlot(allDF,'B2M')

png('figures/dfSOD2.png', width = 1000, height = 1000, res = 300)
VlnPlot(allDF, features = 'SOD2', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'SOD2') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/dfNFKBIA.png', width = 1000, height = 1000, res = 300)
VlnPlot(allDF, features = 'NFKBIA', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'NFKBIA') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/dfB2M.png', width = 1000, height = 1000, res = 300)
VlnPlot(allDF, features = 'B2M', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'B2M') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()