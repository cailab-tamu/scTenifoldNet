SCC <- 'SCC6'
#setwd(paste0('/data/dcosorioh/manuscript/CTX/', SCC))
library(Matrix)
library(scTenifoldNet)
SCC <- 'SCC6'
X <- readMM(paste0(SCC, '_CTL.mtx'))
colnames(X) <- readLines(paste0('barcodes',SCC,'_CTL.txt'))
rownames(X) <- readLines(paste0('genes',SCC,'_CTL.txt'))

Y <- readMM(paste0(SCC, '_CTX.mtx'))
colnames(Y) <- readLines(paste0('barcodes',SCC,'_CTX.txt'))
rownames(Y) <- readLines(paste0('genes',SCC,'_CTX.txt'))

# ## Seurat
# library(Seurat)
# library(harmony)
# X <- CreateSeuratObject(X, project = 'CTL')
# Y <- CreateSeuratObject(Y, project = 'TT')
# ALL <- merge(X,Y)
# ALL <- NormalizeData(ALL)
# ALL <- ScaleData(ALL)
# ALL <- FindVariableFeatures(ALL)
# ALL <- RunPCA(ALL)
# ALL <- RunHarmony(ALL, group.by.vars = 'orig.ident')
# ALL <- RunTSNE(ALL, reduction = 'harmony')
# ALL <- RunUMAP(ALL, reduction = 'harmony', dims = 1:20)

source('https://raw.githubusercontent.com/dosorio/utilities/master/singleCell/scQC.R')
X <- scQC(X)
Y <- scQC(Y)

X <- X[!grepl('^RPL|^RPS|^MT-', rownames(X)),]
Y <- Y[!grepl('^RPL|^RPS|^MT-', rownames(Y)),]

O <- scTenifoldNet(X,Y)
save(O, file = paste0(SCC,'_2.RData'))
