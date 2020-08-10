library(enrichR)
library(scTenifoldNet)
source('https://raw.githubusercontent.com/dosorio/utilities/master/idConvert/hsa2mmu_SYMBOL.R')

nDim <- 30
fileList <- list.files(path = 'results/',pattern = 'sym', full.names = TRUE)

sapply(fileList, function(X){
  dC <- read.csv(X, row.names = 1, stringsAsFactors = FALSE)
  sGenes <- dC$gene[(dC$p.adj < 0.05)]
  sGenes <- toupper(sGenes)
  sGenes <- gsub('^MT-','',sGenes)
  writeLines(sGenes, paste0('geneLists/G_',basename(X)))
  A <- enrichr(sGenes, c('KEGG_2019_Mouse', 'BioPlanet_2019', 'Reactome_2016'))
  A <- lapply(names(A), function(Z){
    DB <- unlist(strsplit(Z,'_'))[1]
    data.frame(DB=DB,A[[Z]])
  })
  A <- do.call(rbind.data.frame, A)
  A <- A[A$Adjusted.P.value < 0.05,]
  A <- A[order(A$Adjusted.P.value),]
  if(! X %in% c('sym10X500DF_Itensor_Dalignment.csv','10X500SCC6_Itensor_Dalignment.csv')){
    A$Genes <- unlist(lapply(A$Genes, function(G){paste0(hsa2mmu_SYMBOL(unlist(strsplit(G,';'))), collapse = ', ')}))
  } else {
    A$Genes <- unlist(lapply(A$Genes, function(G){paste0(unlist(strsplit(G,';')), collapse = ', ')}))
  }
  write.csv(A, file = paste0('results/annotations/A_',basename(X)))
})
