library(Matrix)
library(GGally)
library(fgsea)

KEGG <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Mouse')

load('CCA_morphineNeuron_new.Rdata')
MA <- read.csv('sym10X500morphineNeuron_Itensor_Dalignment.csv', row.names = 1)

zMA <- MA$Z
names(zMA) <- MA$gene

zCCA30 <- dRCCA_d30$Z
names(zCCA30) <- dRCCA_d30$gene

sGenes <- intersect(names(zMA), names(zCCA30))
ALL <- data.frame(MA=zMA[sGenes], CCA30 = zCCA30[sGenes])
rALL <- as.data.frame(apply(ALL,2,rank))

cV <- paste0('rho == ', round(cor(ALL$MA, ALL$CCA30, method = 'sp'),3))
png('MA_CCA.png', width = 1500, height = 1500, res = 300)
ggplot(ALL, aes(MA, CCA30)) + 
  geom_point(color = densCols(ALL, colramp = hcl.colors)) + 
  theme_bw() + 
  xlab(Z-score~(Manifold~Distance)) + 
  ylab(Z-score~(CCA~Distance)) + 
  labs(title = 'CCA vs. MA', subtitle = parse(text = cV)) +
  theme(plot.title = element_text(face = 2))
dev.off()

names(zCCA30) <- toupper(names(zCCA30))
E <- fgseaMultilevel(KEGG,zCCA30, eps = 0)
E$leadingEdge <- unlist(lapply(E$leadingEdge, function(X){paste0(X, collapse = ';')}))
write.csv(E, 'eCCA.csv')       
