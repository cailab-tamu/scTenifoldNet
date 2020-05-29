library(fgsea)
gSets <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
reactomeDB <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
library(ggplot2)
dC <- read.csv('results/sym_10x500SCC6_Itensor_Dalignment.csv')
Z <- dC$Z
names(Z) <- toupper(gsub('MT-','',dC$gene))

E <- fgsea(gSets, stats = Z, nperm = 1e6)
E2 <- fgsea(reactomeDB, stats = Z, nperm = 1e6)

png(filename = 'figures/p1_SCC6.png',width = 1000, height = 1000, pointsize = 30, res = 300)
plotEnrichment(pathway = gSets$`EGFR1 pathway`, stats = Z) + labs(title='EGFR1 pathway', subtitle = paste0('P = ',formatC(E$pval[E$pathway == 'EGFR1 pathway'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p2_SCC6.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = gSets$`Cyclin A-Cdk2-associated events at S phase entry`, stats = Z) + labs(title='Cyclin A-Cdk2-associated events at S phase entry', subtitle = paste0('P = ',formatC(E$pval[E$pathway == 'Cyclin A-Cdk2-associated events at S phase entry'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p3_SCC6.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = reactomeDB$`TP53 Regulates Transcription of Genes Involved in G1 Cell Cycle Arrest Homo sapiens R-HSA-6804116`, stats = Z) + labs(title='G1 Cell Cycle Arrest', subtitle = paste0('P = ',formatC(E2$pval[E2$pathway == 'TP53 Regulates Transcription of Genes Involved in G1 Cell Cycle Arrest Homo sapiens R-HSA-6804116'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p4_SCC6.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = reactomeDB$`Cell Cycle Checkpoints Homo sapiens R-HSA-69620`, stats = Z) + labs(title='Cell Cycle Checkpoints', subtitle = paste0('P = ',formatC(E2$pval[E2$pathway == 'Cell Cycle Checkpoints Homo sapiens R-HSA-69620'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

