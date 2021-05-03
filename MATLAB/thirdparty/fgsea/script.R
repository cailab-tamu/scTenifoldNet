MS <- read.csv('input.txt', row.names = 1)
BC <- MASS::boxcox(MS$drdist~1,plotit=FALSE)
Z <- MS$drdist^abs(BC$x[which.max(BC$y)])
names(Z) <- MS$genelist

library(fgsea)
KEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Human')
BIOP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=BioPlanet_2019')
GOBP <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2018')
REACTOME <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2016')

set.seed(1)
eKEGG <- fgseaMultilevel(KEGG, Z, scoreType = "pos")
set.seed(1)
eBIOP <- fgseaMultilevel(BIOP, Z, scoreType = "pos")
set.seed(1)
eGOBP <- fgseaMultilevel(GOBP, Z, scoreType = "pos")
set.seed(1)
eREACTOME <- fgseaMultilevel(REACTOME, Z, scoreType = "pos")

eKEGG <- eKEGG[eKEGG$ES > 0 & eKEGG$padj < 0.05,]
eBIOP <- eBIOP[eBIOP$ES > 0 & eBIOP$padj < 0.05,]
eGOBP <- eGOBP[eGOBP$ES > 0 & eGOBP$padj < 0.05,]
eREACTOME <- eREACTOME[eREACTOME$ES > 0 & eREACTOME$padj < 0.05,]

E <- do.call(rbind.data.frame, list(eKEGG, eBIOP, eGOBP, eREACTOME))
E$leadingEdge <- unlist(lapply(E$leadingEdge, function(X){paste0(X, collapse = ';')}))
E <- E[order(E$padj, decreasing = FALSE),]
write.csv(E, 'output.txt')
