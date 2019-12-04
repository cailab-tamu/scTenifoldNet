gSets <- read.csv('bioPlanet2019.csv', stringsAsFactors = FALSE)
pIDs <- unique(gSets$PATHWAY_NAME)
gSets <- lapply(pIDs, function(X){
  gSets$GENE_SYMBOL[gSets$PATHWAY_NAME %in% X]
})
names(gSets) <- pIDs

library(fgsea)
library(scTenifoldNet)
library(ggplot2)
mA <- read.csv('results/nonmit_10X500ANEURON_Itensor_Dalignment.csv', row.names = 1)[,1:30]
rownames(mA) <- make.unique(toupper(rownames(mA)))
dC <- dCoexpression(mA, minFC = 0)
Z <- dC$Z
names(Z) <- dC$gene

E <- fgsea(gSets, stats = Z, nperm = 1e6)

png(filename = 'figures/p1_-mtAging.png',width = 1000, height = 1000, pointsize = 30, res = 300)
plotEnrichment(pathway = gSets$`Electron transport chain`, stats = Z) + labs(title='Electron transport\nchain', subtitle = paste0('P = ',formatC(E$pval[E$pathway == 'Electron transport chain'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p2_-mtAging.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = gSets$`Oxidative phosphorylation`, stats = Z) + labs(title='Oxidative\nphosphorylation', subtitle = paste0('P = ',formatC(E$pval[E$pathway == 'Oxidative phosphorylation'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p3_-mtAging.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = gSets$`Parkinson's disease`, stats = Z) + labs(title="Parkinson's disease", subtitle = paste0('P = ',formatC(E$pval[E$pathway == "Parkinson's disease"], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p4_-mtAging.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = gSets$`Nitric oxide signaling pathway`, stats = Z) + labs(title="Nitric oxide\nsignaling pathway", subtitle = paste0('P = ',formatC(E$pval[E$pathway == "Nitric oxide signaling pathway"], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

