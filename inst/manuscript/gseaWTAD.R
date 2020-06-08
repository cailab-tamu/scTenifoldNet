library(fgsea)
library(ggplot2)
load('datasets/WTAD/WTAD.RData')
Z <- WTAD$diffRegulation$Z
names(Z) <- toupper(WTAD$diffRegulation$gene)


gSets <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
reactomeDB <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')

E1 <- fgsea(gSets, stats = Z, nperm = 1e6)
E2 <- fgsea(reactomeDB, stats = Z, nperm = 1e6)


png(filename = 'figures/p1_WTAD.png',width = 1000, height = 1000, pointsize = 30, res = 300)
plotEnrichment(pathway = gSets$`Serotonin HTR1 group and FOS pathway`, stats = Z) + labs(title='Serotonin HTR1 group\nand FOS pathway', subtitle = paste0('P = ',formatC(E1$pval[E1$pathway == 'Serotonin HTR1 group and FOS pathway'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p2_WTAD.png',width = 1000, height = 1000, pointsize = 30, res = 300)
plotEnrichment(pathway = gSets$`Integrin signaling pathway`, stats = Z) + labs(title='Integrin \nsignaling pathway', subtitle = paste0('P = ',formatC(E1$pval[E1$pathway == 'Integrin signaling pathway'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p3_WTAD.png',width = 1000, height = 1000, pointsize = 30, res = 300)
plotEnrichment(pathway = reactomeDB$`Glutamate Neurotransmitter Release Cycle Homo sapiens R-HSA-210500`, stats = Z) + labs(title='Glutamate neurotrans-\nmitter release cycle', subtitle = paste0('P = ',formatC(E1$pval[E1$pathway == 'Glutamate Neurotransmitter Release Cycle Homo sapiens R-HSA-210500'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()
