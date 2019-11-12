library(enrichR)
library(scTenifoldNet)

nDim <- 30
fileList <- list.files(path = 'results/',pattern = 'AN', full.names = TRUE)

sapply(fileList, function(X){
  fileContent <- read.csv(X, row.names = 1, stringsAsFactors = FALSE)
  fileContent <- fileContent[,seq_len(nDim)]
  gList <- rownames(fileContent)
  gList <- gList[grepl('X_', gList)]
  gList <- gsub('X_', '', gList)
  nGenes <- length(gList)
  dC <- dCoexpression(fileContent, nGenes, gList)
  Z <- as.vector(scale(dC$distance))
  P <- pnorm(Z, lower.tail = FALSE)
  table(P < 0.05)
  Q <- p.adjust(P, 'fdr')
  table(Q < 0.05)
  sGenes <- as.vector(dC$gene[Q < 0.05])
  sGenes <- toupper(sGenes)
  sGenes <- gsub('^MT-','',sGenes)
  writeLines(sGenes, paste0('geneLists/G_',basename(X)))
  plot(Z, col=ifelse(Q < 0.05,'red','black'), pch = 16)
  A <- enrichr(sGenes, c('KEGG_2019_Mouse', 'BioPlanet_2019', 'Reactome_2016'))
  A <- do.call(rbind.data.frame, A)
  A <- A[A$P.value < 0.05,]
  A <- A[order(A$Adjusted.P.value),]
  write.csv(A, file = paste0('annotations/A_',basename(X)))
})