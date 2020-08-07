library(fgsea)

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

E1 <- fgseaMultilevel(mmuKEGG, zNeuron)
writeE(E1, 'mmuKEGG_morphineNeurons.csv')

E2 <- fgseaMultilevel(BIOP, zNeuron)
writeE(E2, 'bioplanet_morphineNeurons.csv')

E3 <- fgseaMultilevel(REACTOME, zNeuron)
writeE(E3, 'reactome_morphineNeurons.csv')


#### Dermal Fibroblasts ####
hsaDF <- read.csv('results/sym10X500DF_Itensor_Dalignment.csv', row.names = 1)
hsaDF <- hsaDF[!grepl('^RPL|^RPS|^RP[[:digit:]]', hsaDF$gene, ignore.case = TRUE),]
zDF <- hsaDF$Z
names(zDF) <- toupper(hsaDF$gene)

E1 <- fgseaMultilevel(mmuKEGG, zNeuron)
writeE(E1, 'mmuKEGG_morphineNeurons.csv')

E2 <- fgseaMultilevel(BIOP, zNeuron)
writeE(E2, 'bioplanet_morphineNeurons.csv')

E3 <- fgseaMultilevel(REACTOME, zNeuron)
writeE(E3, 'reactome_morphineNeurons.csv')
