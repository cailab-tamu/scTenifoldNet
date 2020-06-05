library(fgsea)
library(ggplot2)
load('datasets/NKX2-1/Nkx21_AT1.RData')
Z <- O$diffRegulation$Z
names(Z) <- toupper(O$diffRegulation$gene)

CT <- read.csv('datasets/NKX2-1/PanglaoDB_markers_27_Mar_2020.tsv.gz', sep = '\t', stringsAsFactors = FALSE)
CT <- CT[grepl('Mm',CT$species),]
ctNames <- unique(CT$cell.type)
CT <- lapply(ctNames, function(X){
  unique(CT$official.gene.symbol[CT$cell.type %in% X])
})
names(CT) <- ctNames

TISSUE <- read.csv('datasets/NKX2-1/PanglaoDB_markers_27_Mar_2020.tsv.gz', sep = '\t', stringsAsFactors = FALSE)
TISSUE <- TISSUE[grepl('Mm',TISSUE$species),]
tName <- unique(TISSUE$organ)
tName <- tName[!is.na(tName)]
TISSUE <- lapply(tName, function(X){
  unique(TISSUE$official.gene.symbol[TISSUE$organ %in% X])
})
names(TISSUE) <- tName

EC <- fgseaMultilevel(CT, Z)
ET <- fgseaMultilevel(TISSUE, Z)

png(filename = 'figures/p1_NKX21.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(CT$`Pulmonary alveolar type I cells`,Z) + labs(title='Pulmonary alveolar\ntype I cells', subtitle = paste0('P = ',formatC(EC$pval[EC$pathway == 'Pulmonary alveolar type I cells'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p2_NKX21.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(TISSUE[['GI tract']],Z) + labs(title='GI tract', subtitle = paste0('P = ',formatC(ET$pval[ET$pathway == 'GI tract'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

