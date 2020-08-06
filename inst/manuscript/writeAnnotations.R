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

mmuNeuron <- read.csv('results/sym10X500morphineNeuron_Itensor_Dalignment.csv', row.names = 1)
mmuNeuron <- mmuNeuron[!grepl('^RPL|^RPS|^RP[[:digit:]]', mmuNeuron$gene, ignore.case = TRUE),]
zNeuron <- mmuNeuron$Z
names(zNeuron) <- toupper(mmuNeuron$gene)

E1 <- fgseaMultilevel(mmuKEGG, zNeuron)
E1 <- E1[E1$NES > 0 & E1$padj < 0.05,]
E1 <- E1[order(E1$NES, decreasing = TRUE),]

E2 <- fgseaMultilevel(BIOP, zNeuron)
E2 <- E2[E2$NES > 0 & E2$padj < 0.05,]
E2 <- E2[order(E2$NES, decreasing = TRUE),]

E3 <- fgseaMultilevel(REACTOME, zNeuron)
E3 <- E3[E3$NES > 0 & E3$padj < 0.05,]
E3 <- E3[order(E3$NES, decreasing = TRUE),]
