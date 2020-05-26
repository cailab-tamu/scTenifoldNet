library(Matrix)
library(ggplot2)
library(MAST)
library(Seurat)
library(patchwork)
library(ggrepel)

X <- readMM('datasets/morphineNeurons/mock/control.mtx')
colnames(X) <- readLines('datasets/morphineNeurons/mock/controlBarcodes.tsv')
rownames(X) <- readLines('datasets/morphineNeurons/mock/controlGenes.tsv')

Y <- readMM('datasets/morphineNeurons/morphine/morphine.mtx')
colnames(Y) <- readLines('datasets/morphineNeurons/morphine/morphineBarcodes.tsv')
rownames(Y) <- readLines('datasets/morphineNeurons/morphine/morphineGenes.tsv')

avgExpressionX <- rowMeans(X)
avgExpressionY <- rowMeans(Y)

X <- CreateSeuratObject(counts = X)
Y <- CreateSeuratObject(counts = Y)
ALL <- merge(X,Y)
Idents(ALL) <- ifelse(grepl('C',Idents(ALL)),'Mock','Morphine')
ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)
DE <- FindMarkers(ALL, ident.1 = 'Mock', ident.2 = 'Morphine', test.use = 'MAST', logfc.threshold = 0, min.pct = 0.05)

DR <- read.csv('results/sym10X500morphineNeuron_Itensor_Dalignment.csv', row.names = 1)
DE <- DE[DR$gene,]

plotData <- data.frame(X=avgExpressionX[DR$gene], Y= avgExpressionY[DR$gene], DR = DR$FC, DE = DE$avg_logFC, geneName = DR$gene)
plotData$geneName[plotData$geneName %in% DR$gene[!DR$p.adj < 0.05]] <- NA

A <- ggplot(plotData, aes(log1p(X),log1p(DR), label = geneName)) + geom_point(col=densCols(log1p(plotData[,c('X','DR')]))) + ylab(expression(log('Fold Change Distance Mock - Morphine')))+ xlab(expression(log('Mean Gene Expression + 1'))) + theme_bw() + geom_text_repel() + labs(title = 'Mock Sample')
B <- ggplot(plotData, aes(log1p(Y),log1p(DR), label = geneName)) + geom_point(col=densCols(log1p(plotData[,c('Y','DR')]))) + ylab(expression(log('Fold Change Distance Mock - Morphine')))+ xlab(expression(log('Mean Gene Expression + 1'))) + theme_bw() + geom_text_repel()+ labs(title = 'Morphine Sample')
C <- ggplot(plotData, aes(DE,log1p(DR), label = geneName)) + geom_point(col=densCols(log1p(plotData[,c('DE','DR')]))) + ylab(expression(log('Fold Change Distance Mock - Morphine')))+ xlab(expression(log('Fold Change Mock - Morphine'))) + theme_bw() + geom_text_repel()+ labs(title = 'Differential Expression')

png('ctl_Morphine.png', width = 6000, height = 2000, res = 300)
A + B + C
dev.off()


X <- readMM('datasets/dermalFibroblasts/unstimulated/matrix.mtx')
colnames(X) <- readLines('datasets/dermalFibroblasts/unstimulated/barcodes.tsv')
rownames(X) <- read.table('datasets/dermalFibroblasts/unstimulated/genes.tsv', header = FALSE)[,2]

Y <- readMM('datasets/dermalFibroblasts/stimulated/matrix.mtx')
colnames(Y) <- readLines('datasets/dermalFibroblasts/stimulated/barcodes.tsv')
rownames(Y) <- read.table('datasets/dermalFibroblasts/stimulated/genes.tsv', header = FALSE)[,2]

avgExpressionX <- rowMeans(X)
avgExpressionY <- rowMeans(Y)

X <- CreateSeuratObject(counts = X, project = 'C')
Y <- CreateSeuratObject(counts = Y, project = 'T')
ALL <- merge(X,Y)
ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)
DE <- FindMarkers(ALL, ident.1 = 'C', ident.2 = 'T', test.use = 'MAST', logfc.threshold = 0, min.pct = 0.05)

DR <- read.csv('results/sym10X500DF_Itensor_Dalignment.csv', row.names = 1)
DE <- DE[DR$gene,]

plotData <- data.frame(X=avgExpressionX[DR$gene], Y= avgExpressionY[DR$gene], DR = DR$FC, DE = DE$avg_logFC, geneName = DR$gene)
plotData$geneName[plotData$geneName %in% DR$gene[!DR$p.adj < 0.05]] <- NA

A <- ggplot(plotData, aes(log1p(X),log1p(DR), label = geneName)) + geom_point(col=densCols(log1p(plotData[,c('X','DR')]))) + ylab(expression(log('Fold Change Distance Control - dsRNA')))+ xlab(expression(log('Mean Gene Expression + 1'))) + theme_bw() + geom_text_repel() + labs(title = 'Control Sample')
B <- ggplot(plotData, aes(log1p(Y),log1p(DR), label = geneName)) + geom_point(col=densCols(log1p(plotData[,c('Y','DR')]))) + ylab(expression(log('Fold Change Distance Control - dsRNA')))+ xlab(expression(log('Mean Gene Expression + 1'))) + theme_bw() + geom_text_repel()+ labs(title = 'dsRNA Sample')
C <- ggplot(plotData, aes(DE,log1p(DR), label = geneName)) + geom_point(col=densCols(log1p(plotData[,c('DE','DR')]))) + ylab(expression(log('Fold Change Distance Control - dsRNA')))+ xlab(expression(log('Fold Change Control - dsRNA'))) + theme_bw() + geom_text_repel()+ labs(title = 'Differential Expression')

png('ctl_dsRNA.png', width = 6000, height = 2000, res = 300)
A + B + C
dev.off()


X <- readMM('datasets/agingNeurons/Y/yMatrix.mtx')
colnames(X) <- readLines('datasets/agingNeurons/Y/yBarcodes.txt')
rownames(X) <- readLines('datasets/agingNeurons/Y/yGenes.txt')

Y <- readMM('datasets/agingNeurons/O/oMatrix.mtx')
colnames(Y) <- readLines('datasets/agingNeurons/O/oBarcodes.txt')
rownames(Y) <- readLines('datasets/agingNeurons/O/oGenes.txt')

avgExpressionX <- rowMeans(X)
avgExpressionY <- rowMeans(Y)

X <- CreateSeuratObject(counts = X, project = 'Y')
Idents(X) <- 'Y'
Y <- CreateSeuratObject(counts = Y, project = 'O')
Idents(Y) <- 'O'
ALL <- merge(X,Y)
ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)
DE <- FindMarkers(ALL, ident.1 = 'Y', ident.2 = 'O', test.use = 'MAST', logfc.threshold = 0, min.pct = 0.05)

DR <- read.csv('results/sym10X500ANEURON_Itensor_Dalignment.csv', row.names = 1)
DE <- DE[DR$gene,]

plotData <- data.frame(X=avgExpressionX[DR$gene], Y= avgExpressionY[DR$gene], DR = DR$FC, DE = DE$avg_logFC, geneName = DR$gene)
plotData$geneName[plotData$geneName %in% DR$gene[!DR$p.adj < 0.05]] <- NA

plot(plotData$X, plotData$Y)
A <- ggplot(plotData, aes(log1p(X),log1p(DR), label = geneName)) + geom_point(col=densCols(log1p(plotData[,c('X','DR')]))) + ylab(expression(log('Fold Change Distance Young - Old Mice')))+ xlab(expression(log('Mean Gene Expression + 1'))) + theme_bw() + geom_text_repel() + labs(title = 'Young mice Sample')
B <- ggplot(plotData, aes(log1p(Y),log1p(DR), label = geneName)) + geom_point(col=densCols(log1p(plotData[,c('Y','DR')]))) + ylab(expression(log('Fold Change Distance Young - Old Mice')))+ xlab(expression(log('Mean Gene Expression + 1'))) + theme_bw() + geom_text_repel()+ labs(title = 'Old mice Sample')
C <- ggplot(plotData, aes(DE,log1p(DR), label = geneName)) + geom_point(col=densCols(log1p(plotData[,c('DE','DR')]))) + ylab(expression(log('Fold Change Distance Young - Old Mouse')))+ xlab(expression(log('Fold Change Young - Old Mouse'))) + theme_bw() + geom_text_repel()+ labs(title = 'Differential Expression')

png('ctl_Aging.png', width = 6000, height = 2000, res = 300)
A + B + C
dev.off()

