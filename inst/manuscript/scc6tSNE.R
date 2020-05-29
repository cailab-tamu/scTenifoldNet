library(Seurat)
library(Matrix)
library(ggplot2)
library(scTenifoldNet)
library(RColorBrewer)
library(harmony)

SCC <- 'SCC6'
X <- readMM(paste0('datasets/SCC6/',SCC, '_CTL.mtx'))
colnames(X) <- readLines(paste0('datasets/SCC6/barcodes',SCC,'_CTL.txt'))
rownames(X) <- readLines(paste0('datasets/SCC6/genes',SCC,'_CTL.txt'))

Y <- readMM(paste0('datasets/SCC6/',SCC, '_CTX.mtx'))
colnames(Y) <- readLines(paste0('datasets/SCC6/barcodes',SCC,'_CTX.txt'))
rownames(Y) <- readLines(paste0('datasets/SCC6/genes',SCC,'_CTX.txt'))

X <- CreateSeuratObject(X, project = 'Control')
Y <- CreateSeuratObject(Y, project = 'Cetuximab')

allSCC6 <- merge(X, Y)
allSCC6 <- NormalizeData(allSCC6)
allSCC6 <- ScaleData(allSCC6)
allSCC6 <- FindVariableFeatures(allSCC6)
allSCC6 <- RunPCA(allSCC6, verbose = FALSE)
allSCC6 <- RunHarmony(allSCC6, group.by.vars = 'orig.ident')
allSCC6 <- RunTSNE(allSCC6, reduction = 'harmony')
allSCC6 <- RunUMAP(allSCC6, dims = 1:20, reduction = 'harmony')
UMAPPlot(allSCC6)

gID <- Idents(allSCC6)
levels(gID) <- rev(levels(gID))
df <- data.frame(allSCC6@reductions$tsne@cell.embeddings, Group = as.factor(gID))

png('figures/tsneSCC6.png', width = 1000, height = 1000, res = 300, pointsize = 30)
ggplot(df, mapping = aes(x = tSNE_1, y = tSNE_2, color = Group)) + geom_point(pch = 20, cex = 1)  + theme_bw() + xlab('t-SNE 1') + ylab('t-SNE 2') + scale_color_brewer(palette='Set1') + theme(legend.position = "none")
dev.off()

png('figures/scc6DUSP4.png', width = 1000, height = 1000, res = 300)
VlnPlot(allSCC6, features = 'DUSP4', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'DUSP4') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/scc6ITGA3.png', width = 1000, height = 1000, res = 300)
VlnPlot(allSCC6, features = 'ITGA3', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'ITGA3') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()

png('figures/scc6LIF.png', width = 1000, height = 1000, res = 300)
VlnPlot(allSCC6, features = 'LIF', cols = c(rep('#E41A1C',1), rep("#377EB8",1)), pt.size = 0.25)+ theme_bw() + labs(title = 'LIF') + theme(plot.title = element_text(size=20),legend.position = "none") + xlab('Sample')
dev.off()
