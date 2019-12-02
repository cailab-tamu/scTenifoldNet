source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa_ENSEMBL2SYMBOL.R')
reactomeDB <- read.csv('reactomeDB2019.txt', stringsAsFactors = FALSE, sep = '\t', header = FALSE)
reactomeDB <- reactomeDB[reactomeDB$V6 == 'Homo sapiens',]
symbolID <- hsa_ENSEMBL2SYMBOL(reactomeDB$V1)
pathID <- unique(reactomeDB$V4)
reactomeDB <- lapply(pathID, function(X){
  symbolID$SYMBOL[symbolID$GENEID %in% reactomeDB$V1[reactomeDB$V4 %in% X]]
})
names(reactomeDB) <- pathID

gSets <- read.csv('bioPlanet2019.csv', stringsAsFactors = FALSE)
pIDs <- unique(gSets$PATHWAY_NAME)
gSets <- lapply(pIDs, function(X){
  gSets$GENE_SYMBOL[gSets$PATHWAY_NAME %in% X]
})
names(gSets) <- pIDs

library(fgsea)
library(scTenifoldNet)
library(ggplot2)
mA <- read.csv('results/10X500DF_Itensor_Dalignment.csv', row.names = 1)[,1:30]
rownames(mA) <- make.unique(toupper(rownames(mA)))
dC <- dCoexpression(mA, minDist = 0)
Z <- dC$Z
names(Z) <- gsub('MT-','',dC$gene)

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

