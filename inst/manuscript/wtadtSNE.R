library(Seurat)
library(Matrix)
library(ggplot2)
library(scTenifoldNet)
library(RColorBrewer)

# WTAD
WT <- read.csv('datasets/WTAD/Xw2.txt', header = FALSE)
rownames(WT) <- readLines('datasets/WTAD/genelist.txt')
colnames(WT) <- paste0('WT_',seq_len(ncol(WT)))
WT <- as.matrix(WT)
WT <- scQC(WT)

AD <- read.csv('datasets/WTAD/Xa2.txt', header = FALSE)
rownames(AD) <- readLines('datasets/WTAD/genelist.txt')
colnames(AD) <- paste0('AD_',seq_len(ncol(AD)))
AD <- as.matrix(AD)
AD <- scQC(AD)

AD <- CreateSeuratObject(AD)
WT <- CreateSeuratObject(WT)

ALL <- merge(WT, AD)
ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)
ALL <- FindVariableFeatures(ALL)
ALL <- RunPCA(ALL, verbose = FALSE)
ALL <- RunTSNE(ALL, check_duplicates=FALSE)


gID <- Idents(ALL)
df <- data.frame(ALL@reductions$tsne@cell.embeddings, Group = as.factor(gID))

png('figures/tsneWTAD.png', width = 1000, height = 1000, res = 300, pointsize = 30)
ggplot(df, mapping = aes(x = tSNE_1, y = tSNE_2, color = Group)) + geom_point(pch = 20, cex = 1)  + theme_bw() + xlab('t-SNE 1') + ylab('t-SNE 2') + scale_color_brewer(palette='Set1') + theme(legend.position = "none")
dev.off()

VlnPlot(ALL,'Daam1')
VlnPlot(ALL,'Mark2')
VlnPlot(ALL,'Chl1')

png('figures/wtadDAAM1.png', width = 1000, height = 1000, res = 300)
VlnPlot(ALL, features = 'Daam1', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'Daam1') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/wtadMARK2.png', width = 1000, height = 1000, res = 300)
VlnPlot(ALL, features = 'Mark2', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'Mark2') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/wtadCHL1.png', width = 1000, height = 1000, res = 300)
VlnPlot(ALL, features = 'Chl1', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'Chl1') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()
