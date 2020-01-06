gSets <- read.csv('bioPlanet2019.csv', stringsAsFactors = FALSE)
pIDs <- unique(gSets$PATHWAY_NAME)
gSets <- lapply(pIDs, function(X){
  gSets$GENE_SYMBOL[gSets$PATHWAY_NAME %in% X]
})
names(gSets) <- pIDs

library(fgsea)
library(scTenifoldNet)
library(ggplot2)
mA <- read.csv('results/10X500ANEURON_Itensor_Dalignment.csv', row.names = 1)[,1:30]
rownames(mA) <- make.unique(toupper(rownames(mA)))
dC <- dRegulation(mA, minFC = 0)
Z <- dC$distance
lambda <- seq(-2, 2, 1/100)
lambda <- lambda[lambda != 0]
BC <- MASS::boxcox(Z~1,lambda =lambda)
BC <- BC$x[which.max(BC$y)]
Z <- Z ^ BC
if(BC < 0){
  Z <- 1/Z
}
Z <- scale(Z)
names(Z) <- gsub('MT-','',dC$gene)

E <- fgsea(gSets, stats = Z, nperm = 1e6)

png(filename = 'figures/p1_Aging.png',width = 1000, height = 1000, pointsize = 30, res = 300)
plotEnrichment(pathway = gSets$`Electron transport chain`, stats = Z) + labs(title='Electron transport\nchain', subtitle = paste0('P = ',formatC(E$pval[E$pathway == 'Electron transport chain'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p2_Aging.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = gSets$`Oxidative phosphorylation`, stats = Z) + labs(title='Oxidative\nphosphorylation', subtitle = paste0('P = ',formatC(E$pval[E$pathway == 'Oxidative phosphorylation'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p3_Aging.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = gSets$`Parkinson's disease`, stats = Z) + labs(title="Parkinson's disease", subtitle = paste0('P = ',formatC(E$pval[E$pathway == "Parkinson's disease"], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p4_Aging.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = gSets$`Alzheimer's disease`, stats = Z) + labs(title="Alzheimer's disease", subtitle = paste0('P = ',formatC(E$pval[E$pathway == "Alzheimer's disease"], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

