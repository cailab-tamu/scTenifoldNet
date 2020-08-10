library(fgsea)
set.seed(1)

writeE <- function(X, outFile){
  X <- X[X$NES > 0 & X$padj < 0.05,]
  X <- X[order(X$NES, decreasing = TRUE),]
  X$leadingEdge <- unlist(lapply(X$leadingEdge, function(G){paste0(G, collapse = ';')}))
  write.csv(X, paste0('../manuscript/results/annotations/',outFile), row.names = FALSE)
}
mmuKEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')
hsaKEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Human')
BIOP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')


#### Morphine Neurons ####
mmuNeuron <- read.csv('results/sym10X500morphineNeuron_Itensor_Dalignment.csv', row.names = 1)
mmuNeuron <- mmuNeuron[!grepl('^RPL|^RPS|^RP[[:digit:]]', mmuNeuron$gene, ignore.case = TRUE),]
zNeuron <- mmuNeuron$Z
names(zNeuron) <- toupper(mmuNeuron$gene)

E1 <- fgsea(mmuKEGG, zNeuron, nperm = 1e6)
writeE(E1, 'mmuKEGG_morphineNeurons.csv')

E2 <- fgsea(BIOP, zNeuron, nperm = 1e6)
writeE(E2, 'bioplanet_morphineNeurons.csv')

E3 <- fgsea(REACTOME, zNeuron, nperm = 1e6)
writeE(E3, 'reactome_morphineNeurons.csv')


#### Dermal Fibroblasts ####
hsaDF <- read.csv('results/sym10X500DF_Itensor_Dalignment.csv', row.names = 1)
hsaDF <- hsaDF[!grepl('^RPL|^RPS|^RP[[:digit:]]', hsaDF$gene, ignore.case = TRUE),]
zDF <- hsaDF$Z
names(zDF) <- toupper(hsaDF$gene)

E1 <- fgsea(hsaKEGG, zDF, nperm = 1e6)
writeE(E1, 'hsaKEGG_DF.csv')

E2 <- fgsea(BIOP, zDF, nperm = 1e6)
writeE(E2, 'bioplanet_DF.csv')

E3 <- fgsea(REACTOME, zDF, nperm = 1e6)
writeE(E3, 'reactome_DF.csv')

#### NKX21 ####
mmuNKX21 <- read.csv('results/sym10X500NKX21_Itensor_Dalignment.csv', row.names = 1)
mmuKX21 <- mmuNKX21[!grepl('^RPL|^RPS|^RP[[:digit:]]', mmuNKX21$gene, ignore.case = TRUE),]
zNKX21 <- mmuNKX21$Z
names(zNKX21) <- toupper(mmuNKX21$gene)

E1 <- fgsea(mmuKEGG, zNKX21, nperm = 1e6)
writeE(E1, 'mmuKEGG_NKX21.csv')

E2 <- fgsea(BIOP, zNKX21, nperm = 1e6)
writeE(E2, 'bioplanet_NKX21.csv')

E3 <- fgsea(REACTOME, zNKX21, nperm = 1e6)
writeE(E3, 'reactome_NKX21.csv')


#### CETUXIMAB ####
hsaSCC <- read.csv('results/sym10x500SCC6_Itensor_Dalignment.csv', row.names = 1)
hsaSCC <- hsaSCC[!grepl('^RPL|^RPS|^RP[[:digit:]]', hsaSCC$gene, ignore.case = TRUE),]
zSCC <- hsaSCC$Z
names(zSCC) <- toupper(hsaSCC$gene)

E1 <- fgsea(hsaKEGG, zSCC, nperm = 1e6)
writeE(E1, 'hsaKEGG_SCC.csv')

E2 <- fgsea(BIOP, zSCC, nperm = 1e6)
writeE(E2, 'bioplanet_SCC.csv')

E3 <- fgsea(REACTOME, zSCC, nperm = 1e6)
writeE(E3, 'reactome_SCC.csv')

#### AD ####
mmuAD <- read.csv('results/sym10X500WTAD_Itensor_Dalignment.csv', row.names = 1)
mmuAD <- mmuAD[!grepl('^RPL|^RPS|^RP[[:digit:]]', mmuAD$gene, ignore.case = TRUE),]
zAD <- mmuAD$Z
names(zAD) <- toupper(mmuAD$gene)

E1 <- fgsea(hsaKEGG, zAD, nperm = 1e6)
writeE(E1, 'hsaKEGG_AD.csv')

E2 <- fgsea(BIOP, zAD, nperm = 1e6)
writeE(E2, 'bioplanet_AD.csv')

E3 <- fgsea(REACTOME, zAD, nperm = 1e6)
writeE(E3, 'reactome_SCC.csv')