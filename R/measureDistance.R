measureDistance <- function(manifoldOutput, nGenes, geneList, q = 0.05){
  dMetric <- sapply(seq_len(nGenes), function(G){
    X <- manifoldOutput[G,]
    Y <- manifoldOutput[(G+nGenes),]
    I <- rbind(X,Y)
    O <- dist(I)
    O <- as.numeric(O)
    return(O)
  })
  names(dMetric) <- geneList
  return(dMetric)
}
