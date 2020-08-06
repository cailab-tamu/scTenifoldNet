library(Matrix)
library(scTenifoldNet)
library(GGally)
library(UpSetR)
library(patchwork)
library(ggrepel)
library(circlize)
library(ComplexHeatmap)
library(fgsea)
source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa2mmu_SYMBOL.R')
colSCC = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
KEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')

load('DM100.RData')
J100 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n100 <- nrow(O$diffRegulation)
z100 <- O$diffRegulation$Z
names(z100) <- toupper(O$diffRegulation$gene)
e100 <- fgseaMultilevel(KEGG, z100)

load('DM200.RData')
J200 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n200 <- nrow(O$diffRegulation)
z200 <- O$diffRegulation$Z
names(z200) <- toupper(O$diffRegulation$gene)
e200 <- fgseaMultilevel(KEGG, z200)

load('DM300.RData')
J300 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n300 <- nrow(O$diffRegulation)
z300 <- O$diffRegulation$Z
names(z300) <- toupper(O$diffRegulation$gene)
e300 <- fgseaMultilevel(KEGG, z300)

load('DM400.RData')
J400 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n400 <- nrow(O$diffRegulation)
z400 <- O$diffRegulation$Z
names(z400) <- toupper(O$diffRegulation$gene)
e400 <- fgseaMultilevel(KEGG, z400)

load('DM500.RData')
J500 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n500 <- nrow(O$diffRegulation)
z500 <- O$diffRegulation$Z
names(z500) <- toupper(O$diffRegulation$gene)
e500 <- fgseaMultilevel(KEGG, z500)


load('DM600.RData')
J600 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n600 <- nrow(O$diffRegulation)
z600 <- O$diffRegulation$Z
names(z600) <- toupper(O$diffRegulation$gene)
e600 <- fgseaMultilevel(KEGG, z600)

load('DM700.RData')
J700 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n700 <- nrow(O$diffRegulation)
z700 <- O$diffRegulation$Z
names(z700) <- toupper(O$diffRegulation$gene)
e700 <- fgseaMultilevel(KEGG, z700)

load('DM800.RData')
J800 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n800 <- nrow(O$diffRegulation)
z800 <- O$diffRegulation$Z
names(z800) <- toupper(O$diffRegulation$gene)
e800 <- fgseaMultilevel(KEGG, z800)

load('DM900.RData')
J900 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n900 <- nrow(O$diffRegulation)
z900 <- O$diffRegulation$Z
names(z900) <- toupper(O$diffRegulation$gene)
e900 <- fgseaMultilevel(KEGG, z900)

load('DM1000.RData')
J1000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n1000 <- nrow(O$diffRegulation)
z1000 <- O$diffRegulation$Z
names(z1000) <- toupper(O$diffRegulation$gene)

load('DM2000.RData')
J2000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n2000 <- nrow(O$diffRegulation)
z2000 <- O$diffRegulation$Z
names(z2000) <- toupper(O$diffRegulation$gene)

load('DM3000.RData')
J3000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n3000 <- nrow(O$diffRegulation)
z3000 <- O$diffRegulation$Z
names(z3000) <- toupper(O$diffRegulation$gene)

load('DM4000.RData')
J4000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n4000 <- nrow(O$diffRegulation)
z4000 <- O$diffRegulation$Z
names(z4000) <- toupper(O$diffRegulation$gene)

load('DM5000.RData')
J5000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n5000 <- nrow(O$diffRegulation)
z5000 <- O$diffRegulation$Z
names(z5000) <- toupper(O$diffRegulation$gene)

load('DM6000.RData')
J6000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n6000 <- nrow(O$diffRegulation)
z6000 <- O$diffRegulation$Z
names(z6000) <- toupper(O$diffRegulation$gene)

load('DM7000.RData')
J7000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000 <- nrow(O$diffRegulation)
z7000 <- O$diffRegulation$Z
names(z7000) <- toupper(O$diffRegulation$gene)
e7000 <- fgseaMultilevel(KEGG, z7000)

load('DM7000x200.RData')
J7000X200 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000X200 <- nrow(O$diffRegulation)
z7000X200 <- O$diffRegulation$Z
names(z7000X200) <- toupper(O$diffRegulation$gene)
e7000x200 <- fgseaMultilevel(KEGG, z7000X200)

load('DM7000x400.RData')
J7000X400 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000X400 <- nrow(O$diffRegulation)
z7000X400 <- O$diffRegulation$Z
names(z7000X400) <- toupper(O$diffRegulation$gene)

load('DM7000x500.RData')
J7000X500 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000X500 <- nrow(O$diffRegulation)
z7000X500 <- O$diffRegulation$Z
names(z7000X500) <- toupper(O$diffRegulation$gene)

load('DM7000x600.RData')
J7000X600 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000X600 <- nrow(O$diffRegulation)
z7000X600 <- O$diffRegulation$Z
names(z7000X600) <- toupper(O$diffRegulation$gene)

load('DM7000x700.RData')
J7000X700 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000X700 <- nrow(O$diffRegulation)
z7000X700 <- O$diffRegulation$Z
names(z7000X700) <- toupper(O$diffRegulation$gene)

load('DM7000x800.RData')
J7000X800 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000X800 <- nrow(O$diffRegulation)
z7000X800 <- O$diffRegulation$Z
names(z7000X800) <- toupper(O$diffRegulation$gene)

load('DM7000x900.RData')
J7000X900 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000X900 <- nrow(O$diffRegulation)
z7000X900 <- O$diffRegulation$Z
names(z7000X900) <- toupper(O$diffRegulation$gene)

load('DM7000x1000.RData')
J7000X1000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000X1000 <- nrow(O$diffRegulation)
z7000X1000 <- O$diffRegulation$Z
names(z7000X1000) <- toupper(O$diffRegulation$gene)

load('DM7000x2000.RData')
J7000X2000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000X2000 <- nrow(O$diffRegulation)
z7000X2000 <- O$diffRegulation$Z
names(z7000X2000) <- toupper(O$diffRegulation$gene)

load('DM7000x3000.RData')
J7000X3000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000X3000 <- nrow(O$diffRegulation)
z7000X3000 <- O$diffRegulation$Z
names(z7000X3000) <- toupper(O$diffRegulation$gene)

load('DM7000x4000.RData')
J7000X4000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000X4000 <- nrow(O$diffRegulation)
z7000X4000 <- O$diffRegulation$Z
names(z7000X4000) <- toupper(O$diffRegulation$gene)

load('DM7000x5000.RData')
J7000X5000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
n7000X5000 <- nrow(O$diffRegulation)
z7000X5000 <- O$diffRegulation$Z
names(z7000X5000) <- toupper(O$diffRegulation$gene)

# load('DM7000x6000.RData')
# J7000X6000 <- O$diffRegulation$gene[O$diffRegulation$p.adj < 0.05]
# n7000X6000 <- nrow(O$diffRegulation)
# z7000X6000 <- O$diffRegulation$Z
# names(z7000X6000) <- O$diffRegulation$gene


sGenes <- intersect(names(z100), names(z200))
sGenes <- intersect(sGenes, names(z300))
sGenes <- intersect(sGenes, names(z400))
sGenes <- intersect(sGenes, names(z500))
sGenes <- intersect(sGenes, names(z600))
sGenes <- intersect(sGenes, names(z700))
sGenes <- intersect(sGenes, names(z800))
sGenes <- intersect(sGenes, names(z900))
sGenes <- intersect(sGenes, names(z1000))
sGenes <- intersect(sGenes, names(z2000))
sGenes <- intersect(sGenes, names(z3000))
sGenes <- intersect(sGenes, names(z4000))
sGenes <- intersect(sGenes, names(z5000))
sGenes <- intersect(sGenes, names(z6000))
sGenes <- intersect(sGenes, names(z7000))
sGenes <- intersect(sGenes, names(z7000X200))
sGenes <- intersect(sGenes, names(z7000X400))
sGenes <- intersect(sGenes, names(z7000X500))
sGenes <- intersect(sGenes, names(z7000X600))
sGenes <- intersect(sGenes, names(z7000X700))
sGenes <- intersect(sGenes, names(z7000X800))
sGenes <- intersect(sGenes, names(z7000X900))
sGenes <- intersect(sGenes, names(z7000X1000))
sGenes <- intersect(sGenes, names(z7000X2000))
sGenes <- intersect(sGenes, names(z7000X3000))
sGenes <- intersect(sGenes, names(z7000X4000))
sGenes <- intersect(sGenes, names(z7000X5000))


allZ <- data.frame(c100x100 = z100[sGenes], c200x200 = z200[sGenes], c300x300 = z300[sGenes], c400x400 = z400[sGenes], c500x500 = z500[sGenes], 
                   c600x600 = z600[sGenes], c700x700 = z700[sGenes], c800x800 = z800[sGenes], c900x900 = z900[sGenes],
                   c1000x1000 = z1000[sGenes], c2000x2000 = z2000[sGenes], c3000x3000 = z3000[sGenes], c4000x4000 = z4000[sGenes], 
                   c5000x5000 = z5000[sGenes], c6000X6000 = z6000[sGenes], c7000x7000 = z7000[sGenes], c7000x200 = z7000X200[sGenes],
                   c7000x400 = z7000X400[sGenes], c7000x500 = z7000X500[sGenes], c7000x600 = z7000X600[sGenes], c7000x800 = z7000X800[sGenes], 
                   c7000x1000 = z7000X1000[sGenes], c7000x2000 = z7000X2000[sGenes], c7000x3000 = z7000X3000[sGenes],
                   c7000x4000 = z7000X4000[sGenes], c7000x5000 = z7000X5000[sGenes])

#rZ <- apply(allZ,2,rank)
PCA <- prcomp(t(allZ))$x
PCA <- as.data.frame(PCA)
PCA$Cells <- gsub('c','',rownames(PCA))
A <- ggplot(PCA, aes(PC1,PC2, label=Cells)) + geom_point() + geom_text_repel(segment.alpha = 0.4, segment.size = .5) + theme_bw()
A

colnames(allZ) <- gsub('c','', colnames(allZ))
rZ <- apply(allZ,2,rank)
allZ <- allZ[order(rowMeans(rZ), decreasing = TRUE),]
MA <- KEGG$`Morphine addiction`[KEGG$`Morphine addiction` %in% rownames(allZ)]
png('rankingHM.png', width = 1500, height = 3500, res = 300)
Heatmap(as.matrix(allZ), show_row_dend = FALSE, show_row_names = FALSE, show_column_dend = FALSE, column_order = seq_len(ncol(allZ)), row_order = seq_len(nrow(allZ)), name = 'Z') +
  rowAnnotation(link = anno_mark(at = which(rownames(allZ) %in% MA), labels = hsa2mmu_SYMBOL(rownames(allZ)[rownames(allZ) %in% MA])))
dev.off()

cValues <- cor(allZ, method = 'sp')
#diag(cValues) <- NA
#cValues <- cValues/max(abs(cValues), na.rm = TRUE)
#rownames(cValues) <- colnames(cValues) <- gsub('c','',colnames(allZ))
#ComplexHeatmap::Heatmap(cValues, name = 'SCC', col = colSCC) #row_order = seq_len(ncol(allZ)), column_order = seq_len(ncol(allZ))
cValues <- reshape2::melt(cValues)
cValues$Var1 <- as.factor(cValues$Var1)
cValues$Var2 <- as.factor(cValues$Var2)
colnames(cValues)[3] <- 'SCC'
B <- ggplot(cValues, aes(Var1, Var2, fill=SCC)) + geom_tile() + theme_minimal() + xlab('Number of Cells') + ylab('Number of Cells') +
  scale_fill_gradientn(colours = c("blue", "white", "red"), values = c(-1,0,1), limits=c(0,1)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


allJ <- fromList(list(c100 = J100, c200 = J200, c300 = J300, c400 = J400, c500 = J500,
                      c600 = J600, c700 = J700, c800 = J800, c900 = J900, c1000 = J1000,
                      c2000 = J2000, c3000 = J3000, c4000 = J4000, c5000 = J5000, c6000 = J6000, c7000 = J7000))
jNames <- gsub('c','',colnames(allJ))
jValues <- sapply(seq_len(ncol(allJ)), function(i){
  sapply(seq_len(ncol(allJ)), function(j){
    tempX <- cbind(allJ[,i], allJ[,j])
    tempX <- rowSums(tempX)
    sum(tempX == 2)/sum(tempX > 0)
  })
})
#jValues[lower.tri(jValues, diag = TRUE)] <- NA
#diag(jValues) <- NA
colnames(jValues) <- rownames(jValues) <- jNames
jValues <- reshape2::melt(jValues)
jValues$Var1 <- as.factor(jValues$Var1)
jValues$Var2 <- as.factor(jValues$Var2)
colnames(jValues)[3] <- 'Jaccard'
C <- ggplot(jValues, aes(Var1, Var2, fill=Jaccard)) + geom_tile() + theme_minimal() + xlab('Number of Cells') + ylab('Number of Cells') + 
  scale_fill_gradientn(colours = c("blue", "white", "red"), values = c(0,0.5,1), limits=c(0,1))


nGenes <- rbind(n100, n200, n300, n400, n500, n600, n700, n800, n900, n1000, n2000, n3000, n4000, n5000, n6000, n7000)
nGenes <- as.data.frame(nGenes)
nGenes$Cells <- gsub('n','',rownames(nGenes))
nGenes$Cells <- factor(nGenes$Cells, levels = nGenes$Cells)
D <- ggplot(nGenes, aes(y = Cells, x = V1)) + 
  geom_bar(stat = 'identity') + 
  theme_bw() + 
  xlab('Number of Tested Genes') + 
  ylab('Number of Cells') + 
  labs(title = 'Mouse Neurons') + 
  theme(plot.title = element_text(face = 2))

png('DownSampling.png', width = 1500, height = 3500, res = 300)
A / B / D
dev.off()
