library(fgsea)
library(Matrix)
library(ggplot2)
library(patchwork)

DR <- read.csv('../manuscript/results/sym10X500morphineNeuron_Itensor_Dalignment.csv', row.names = 1)
Z <- DR$Z
names(Z) <- toupper(DR$gene)

KEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')

E <- fgseaMultilevel(KEGG, Z)
E <- E[E$NES > 0 & E$padj < 0.05,]
E <- E[order(E$padj, decreasing = TRUE),]
E$pathway <- factor(E$pathway, levels = E$pathway)
p0 <- ggplot(E, aes(-log10(padj),pathway, fill = NES)) + geom_bar(stat = 'identity') + ylab('KEGG Pathway') + xlab(expression(-log[10]~Adjusted~P-Value)) + theme_bw()

E <- E[order(E$NES, decreasing = TRUE),]
A <- E[grepl('addic',E$pathway),]
A$pathway <- factor(A$pathway, levels = rev(A$pathway))


pA <- ggplot(A, aes(-log10(padj),pathway)) + geom_bar(stat = 'identity') + ylab('KEGG Pathway') + xlab(expression(-log[10]~Adjusted~P-Value)) + theme_bw()
pB <- ggplot(A, aes(NES,pathway)) + geom_bar(stat = 'identity') + ylab('KEGG Pathway') + xlab('Normalized Enrichment Score') + theme_bw()

png('eAddiction.png', width = 1000, height = 1500, res = 300)
(pA / pB)
dev.off()

png('eAll.png', width = 2000, height = 2500, res = 300)
p0
dev.off()
