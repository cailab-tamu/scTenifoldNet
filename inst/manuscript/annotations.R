library(enrichR)
library(scTenifoldNet)

nDim <- 30
fileList <- list.files(path = 'results/',pattern = 'sym', full.names = TRUE)

sapply(fileList, function(X){
  dC <- read.csv(X, row.names = 1, stringsAsFactors = FALSE)
  sGenes <- dC$gene[(dC$p.adj < 0.05)]
  sGenes <- toupper(sGenes)
  sGenes <- gsub('^MT-','',sGenes)
  writeLines(sGenes, paste0('geneLists/G_',basename(X)))
  plot(dC$Z, col=ifelse(dC$p.adj < 0.05,'red','black'), pch = 16)
  plot(-log10(seq(0,1,length.out = nrow(dC))), -log10(dC$p.value))
  A <- enrichr(sGenes, c('KEGG_2019_Mouse', 'BioPlanet_2019', 'Reactome_2016'))
  A <- do.call(rbind.data.frame, A)
  A <- A[A$Adjusted.P.value < 0.05,]
  A <- A[order(A$Adjusted.P.value),]
  write.csv(A, file = paste0('results/annotations/A_',basename(X)))
})
