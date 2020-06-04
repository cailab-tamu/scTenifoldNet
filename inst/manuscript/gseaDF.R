library(fgsea)
gSets <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
reactomeDB <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
library(ggplot2)
load('datasets/dermalFibroblasts/DF_2.RData')
dC <- DF$diffRegulation
Z <- dC$Z
names(Z) <- toupper(gsub('MT-','',dC$gene))

E <- fgsea(gSets, stats = Z, nperm = 1e6)
E2 <- fgsea(reactomeDB, stats = Z, nperm = 1e6)

png(filename = 'figures/p1_DF.png',width = 1000, height = 1000, pointsize = 30, res = 300)
plotEnrichment(pathway = gSets$`Interferon alpha/beta signaling`, stats = Z) + labs(title='Interferon alpha/beta\nsignaling', subtitle = paste0('P = ',formatC(E$pval[E$pathway == 'Interferon alpha/beta signaling'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p2_DF.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = gSets$`Translation`, stats = Z) + labs(title='Translation', subtitle = paste0('P = ',formatC(E$pval[E$pathway == 'Translation'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()
