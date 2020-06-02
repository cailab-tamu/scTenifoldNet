library(Matrix)
library(Seurat)

WT <- Matrix(as.matrix(read.csv('Xw2.txt', header = FALSE)))
AD <- Matrix(as.matrix(read.csv('Xa2.txt', header = FALSE)))
rownames(WT) <- rownames(AD) <- readLines('genelist.txt')
colnames(WT) <- paste0('C',seq_len(ncol(WT)))
colnames(AD) <- paste0('C',seq_len(ncol(AD)))

WT <- CreateSeuratObject(WT, project = 'WT')
AD <- CreateSeuratObject(AD, project = 'AD')

ALL <- merge(WT,AD)
ALL <- NormalizeData(ALL)
ALL <- ScaleData(ALL)
ALL <- FindVariableFeatures(ALL)
ALL <- RunPCA(ALL)
ALL <- RunTSNE(ALL, check_duplicates = FALSE)
TSNEPlot(ALL)

source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/plotDR.R')
DR <- read.csv('results_WTAD.csv', row.names = 1, stringsAsFactors = FALSE)
plotDR(list(diffRegulation=DR))

Z <- DR$FC
names(Z) <- toupper(DR$gene)
library(fgsea)
KEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
BIOP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
E <- fgseaMultilevel(KEGG,Z)
E <- E[E$padj < 0.05,]
