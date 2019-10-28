#' @export dCoexpression
#' @importFrom stats dist
#' @title dCoexpression
#' @description ...
#' @param manifoldOutput ...
#' @param nGenes ...
#' @param geneList ...
#' @param q ...
#' @return ...
#' @references ...

# manifoldOutput <- read.csv('Z:/Cailab/aging_mouse_brain/microglia/Y_O_ManifoldAlignment_microglia.csv', row.names = 1)
# a0 <- read.csv('Z:/Cailab/mouse_morphine_glia_GSE118918_RAW/aln0.txt', header = FALSE)
# a1 <- read.csv('Z:/Cailab/mouse_morphine_glia_GSE118918_RAW/aln1.txt', header = FALSE)
# geneList <- rownames(manifoldOutput)
# geneList <- geneList[grepl('Y', geneList)]
# geneList <- unlist(lapply(strsplit(geneList, '_'), function(X){X[2]}))
# nGenes <- length(geneList)
#manifoldOutput <- rbind(a0,a1)

dCoexpression <- function(manifoldOutput, nGenes, geneList){
  dMetric <- sapply(seq_len(nGenes), function(G){
    X <- manifoldOutput[G,]
    Y <- manifoldOutput[(G+nGenes),]
    I <- rbind(X,Y)
    O <- dist(I)
    O <- as.numeric(O)
    return(O)
  })
  
  dMetric[is.na(dMetric)] <- 0
  dMetric[dMetric == 0] <- 1e-16
  bcExp <- MASS::boxcox(dMetric~1, plotit=FALSE)
  bcExp <- abs(bcExp$x[which.max(bcExp$y)])
  
  if(bcExp != 0){
    nDist <- dMetric ^ bcExp  
  } else {
    nDist <- dMetric
  }
  
  sDist <- scale(nDist)
  
  pValues <- pnorm(sDist, lower.tail = FALSE)
  pAdjusted <- p.adjust(pValues, method = 'fdr')
  dOut <- data.frame(
    gene = geneList, 
    distance = dMetric,
    Z = sDist,
    p.value = pValues,
    p.adj = pAdjusted
  )
  dOut <- dOut[order(dOut$p.value),]
  return(dOut)
}

# dCoex <- dCoexpression(manifoldOutput = manifoldOutput, nGenes = nGenes, geneList = geneList)
# dCoex <- dCoex[dCoex$Z >= 0,]
# 
# library(fgsea)
# C2 <- gmtPathways('~/../Downloads/c2.all.v7.0.symbols.gmt')
# Z <- dCoex$Z
# names(Z) <- toupper(dCoex$gene)
# 
# o <- fgsea(C2, stats = Z, nperm = 1e3)
# o <- o[o$NES > 0,]
# o <- o[order(o$NES, decreasing = TRUE),]
