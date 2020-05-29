library(Matrix)
setwd('/data/dcosorioh/manuscript/NKX2-1/Data')

# WT <- readMM('GSM3716703_Nkx2-1_control_scRNAseq_matrix.mtx.gz')
# colnames(WT) <- readLines('GSM3716703_Nkx2-1_control_scRNAseq_barcodes.tsv.gz')
# rownames(WT) <- read.table('GSM3716703_Nkx2-1_control_scRNAseq_genes.tsv.gz', sep = '\t', header = FALSE, stringsAsFactors = FALSE)[,2]
# KO <- readMM('GSM3716704_Nkx2-1_mutant_scRNAseq_matrix.mtx.gz')
# colnames(KO) <- readLines('GSM3716704_Nkx2-1_mutant_scRNAseq_barcodes.tsv.gz')
# rownames(KO) <- read.table('GSM3716704_Nkx2-1_mutant_scRNAseq_genes.tsv.gz', sep = '\t', header = FALSE, stringsAsFactors = FALSE)[,2]

# source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')
# WT <- scQC(WT)
# KO <- scQC(KO)
# 
# WT <- WT[!grepl('^Rpl|^Rps|^Mt-',rownames(WT), ignore.case = TRUE),]
# KO <- KO[!grepl('^Rpl|^Rps|^Mt-',rownames(KO), ignore.case = TRUE),]

## Seurat
# library(Seurat)
# library(harmony)
# WT <- CreateSeuratObject(WT, project = 'WT')
# KO <- CreateSeuratObject(KO, project = 'KO')
# ALL <- merge(WT,KO)
# ALL <- NormalizeData(ALL)
# ALL <- ScaleData(ALL)
# ALL <- FindVariableFeatures(ALL)
# ALL <- RunPCA(ALL)
# ALL <- RunHarmony(ALL, group.by.vars = 'orig.ident')
# ALL <- RunTSNE(ALL, reduction = 'harmony')
# ALL <- RunUMAP(ALL, reduction = 'harmony', dims = 1:20)
# png('NKX21.png', width = 2000, height = 1500, res = 300)
# TSNEPlot(ALL)
# dev.off()
# 
# ALL <- FindNeighbors(ALL, reduction = 'harmony')
# ALL <- FindClusters(ALL, resolution = 0.05)
# TSNEPlot(ALL)
# 
# #DE <- FindAllMarkers(ALL, test.use = 'MAST', logfc.threshold = .1)
# #write.csv(DE, file = 'deNKX21.csv')
# DE <- read.csv('deNKX21.csv', row.names = 1, stringsAsFactors = FALSE)
# library(fgsea)
# CT <- read.csv('~/../Downloads/PanglaoDB_markers_27_Mar_2020.tsv.gz', sep = '\t', stringsAsFactors = FALSE)
# CT <- CT[CT$organ %in% 'Lungs',]
# CTnames <- unique(CT$cell.type)
# CT <- lapply(CTnames, function(X){
#   unique(CT$official.gene.symbol[CT$cell.type %in% X])
# })
# names(CT) <- CTnames
# dCT <- sapply(seq_along(unique(Idents(ALL)))-1, function(X){
#   tDE <- DE[DE$cluster == X,]
#   Z <- tDE$avg_logFC
#   names(Z) <- toupper(tDE$gene)
#   tCT <- fgseaMultilevel(CT, Z)
#   tCT <- tCT[order(tCT$padj),]
#   tCT <- tCT[tCT$NES > 0 & tCT$padj < 0.05,]
#   tCT$pathway[1]
# })
# levels(Idents(ALL)) <- dCT
# TSNEPlot(ALL)
# 
# table(Idents(ALL), ALL$orig.ident)
# 
# wtAT1 <- ALL@assays$RNA@counts[,Idents(ALL) %in% 'Pulmonary alveolar type I cells' & ALL$orig.ident %in% 'WT']
# writeMM(wtAT1, file = 'wtAT1.mtx')
# writeLines(rownames(wtAT1), 'genes_wtAT1.txt')
# writeLines(colnames(wtAT1), 'barcodes_wtAT1.txt')
# koAT1 <- ALL@assays$RNA@counts[,Idents(ALL) %in% 'Pulmonary alveolar type I cells' & ALL$orig.ident %in% 'KO']
# writeMM(koAT1, file = 'koAT1.mtx')
# writeLines(rownames(koAT1), 'genes_koAT1.txt')
# writeLines(colnames(koAT1), 'barcodes_koAT1.txt')
# wtAT2 <- ALL@assays$RNA@counts[,Idents(ALL) %in% 'Pulmonary alveolar type II cells' & ALL$orig.ident %in% 'WT']
# writeMM(wtAT2, file = 'wtAT2.mtx')
# writeLines(rownames(wtAT2), 'genes_wtAT2.txt')
# writeLines(colnames(wtAT2), 'barcodes_wtAT2.txt')
# koAT2 <- ALL@assays$RNA@counts[,Idents(ALL) %in% 'Pulmonary alveolar type II cells' & ALL$orig.ident %in% 'KO']
# writeMM(koAT2, file = 'koAT2.mtx')
# writeLines(rownames(koAT2), 'genes_koAT2.txt')
# writeLines(colnames(koAT2), 'barcodes_koAT2.txt')

library(Matrix)
WT <- readMM('wtAT1.mtx')
rownames(WT) <- readLines('genes_wtAT1.txt')
colnames(WT) <- readLines('barcodes_wtAT1.txt')

KO <- readMM('koAT1.mtx')
rownames(KO) <- readLines('genes_koAT1.txt')
colnames(KO) <- readLines('barcodes_koAT1.txt')

source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')
WT <- scQC(WT)
KO <- scQC(KO)

WT <- WT[!grepl('^Rpl|^Rps|^Mt-',rownames(WT), ignore.case = TRUE),]
KO <- KO[!grepl('^Rpl|^Rps|^Mt-',rownames(KO), ignore.case = TRUE),]

library(scTenifoldNet)
O <- scTenifoldNet(WT, KO, )
save(O, file = 'Nkx21_AT1.RData')


