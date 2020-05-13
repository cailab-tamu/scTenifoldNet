library(fgsea)
gSets <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
reactomeDB <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')
library(ggplot2)
dC <- read.csv('results/sym10X500DF_Itensor_Dalignment.csv')
Z <- dC$Z
names(Z) <- toupper(gsub('MT-','',dC$gene))

E <- fgsea(gSets, stats = Z, nperm = 1e6)
E2 <- fgsea(reactomeDB, stats = Z, nperm = 1e6)

png(filename = 'figures/p1_DF.png',width = 1000, height = 1000, pointsize = 30, res = 300)
plotEnrichment(pathway = gSets$`Interferon alpha/beta signaling`, stats = Z) + labs(title='Interferon alpha/beta\nsignaling', subtitle = paste0('P = ',formatC(E$pval[E$pathway == 'Interferon alpha/beta signaling'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p2_DF.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = gSets$`Influenza viral RNA transcription and replication`, stats = Z) + labs(title='Viral RNA transcription\nand replication', subtitle = paste0('P = ',formatC(E$pval[E$pathway == 'Influenza viral RNA transcription and replication'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p3_DF.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = reactomeDB$`Peptide chain elongation`, stats = Z) + labs(title='Peptide chain\nelongation', subtitle = paste0('P = ',formatC(E2$pval[E2$pathway == 'Peptide chain elongation'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

png(filename = 'figures/p4_DF.png',width = 1000, height = 1000, pointsize = 20, res = 300)
plotEnrichment(pathway = reactomeDB$`Infectious disease`, stats = Z) + labs(title='Infectious\ndisease', subtitle = paste0('P = ',formatC(E2$pval[E2$pathway == 'Infectious disease'], digits = 3))) + ylab('Enrichment Score') + xlab('Rank') + theme_bw() + theme(plot.title = element_text(size=20))
dev.off()

