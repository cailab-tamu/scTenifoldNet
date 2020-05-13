library(fgsea)
gSets <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')

library(ggplot2)
dC <- read.csv('results/sym10X500ANEURON_Itensor_Dalignment.csv')
Z <- dC$Z
names(Z) <- toupper(gsub('MT-','',dC$gene))

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

