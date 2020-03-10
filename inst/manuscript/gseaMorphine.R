gSets <- read.csv('bioPlanet2019.csv', stringsAsFactors = FALSE)
pIDs <- unique(gSets$PATHWAY_NAME)
gSets <- lapply(pIDs, function(X){
  gSets$GENE_SYMBOL[gSets$PATHWAY_NAME %in% X]
})
names(gSets) <- pIDs


KEGG <- readLines('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
KEGG <- strsplit(KEGG, '\t')
KEGGnames <- lapply(KEGG, function(X){X[1]})
KEGG <- lapply(KEGG, function(X){X[-c(1,2)]})
names(KEGG) <- unlist(KEGGnames)

library(fgsea)
library(scTenifoldNet)
library(ggplot2)
mA <- read.csv('results/10X500morphineNeuron_Itensor_Dalignment.csv', row.names = 1)[,1:30]
rownames(mA) <- make.unique(toupper(rownames(mA)))
dC <- dRegulation(mA, minFC = 0)
Z <- dC$Z
names(Z) <- dC$gene

E <- fgsea(gSets, stats = Z, nperm = 1e6)

png(filename = 'figures/p1_Morphine.png',width = 1000, height = 1000, pointsize = 30, res = 300)
plotEnrichment(pathway = gSets$`Opioid signaling`, stats = Z) + labs(title='Opioid\nsignaling', subtitle = paste0('P = ',formatC(E$pval[E$pathway == 'Opioid signaling'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p2_Morphine.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = gSets$`G-protein signaling pathways`, stats = Z) + labs(title='G-protein\nsignaling pathways', subtitle = paste0('P = ',formatC(E$pval[E$pathway == 'G-protein signaling pathways'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p3_Morphine.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = gSets$`GABA A and B receptor activation`, stats = Z) + labs(title='GABA A and B\nreceptor activation', subtitle = paste0('P = ',formatC(E$pval[E$pathway == 'GABA A and B receptor activation'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p4_Morphine.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = gSets$`Calcium signaling pathway`, stats = Z) + labs(title='Calcium signaling\npathway', subtitle = paste0('P = ',formatC(E$padj[E$pathway == 'Calcium signaling pathway'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

E <- fgsea(KEGG, stats = Z, nperm = 1e6)
png(filename = 'figures/p5_Morphine.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = KEGG$`Morphine addiction`, stats = Z) + labs(title='Morphine\naddiction', subtitle = paste0('P = ',formatC(E$padj[E$pathway == 'Morphine addiction'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

