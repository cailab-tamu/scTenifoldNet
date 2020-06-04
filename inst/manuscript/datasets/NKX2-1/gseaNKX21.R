library(fgsea)
load('Nkx21_AT1.RData')
Z <- O$diffRegulation$Z
names(Z) <- toupper(O$diffRegulation$gene)

CT <- read.csv('PanglaoDB_markers_27_Mar_2020.tsv.gz', sep = '\t', stringsAsFactors = FALSE)
CT <- CT[grepl('Mm',CT$species),]
ctNames <- unique(CT$cell.type)
CT <- lapply(ctNames, function(X){
  unique(CT$official.gene.symbol[CT$cell.type %in% X])
})
names(CT) <- ctNames

E <- fgseaMultilevel(CT, Z)
plotEnrichment(CT$`Pulmonary alveolar type I cells`,Z)
