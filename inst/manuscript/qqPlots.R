library(Matrix)
library(apcluster)
library(scTenifoldNet)
library(ggplot2)
library(ggrepel)

# load('datasets/dermalFibroblasts/DF_2.RData')
# write.csv(DF$diffRegulation, file = 'results/sym10X500DF_Itensor_Dalignment.csv')

# load('datasets/NKX2-1/Nkx21_AT1.RData')
# write.csv(O$diffRegulation, file = 'results/sym10X500NKX21_Itensor_Dalignment.csv')

fileList <- list.files('results/', pattern = 'sym', full.names = TRUE)
sID <- c('SCC6','Aging', 'DermalFibroblasts', 'Morphine', 'NKX21')

sapply(seq_along(fileList), function(X){
  dC <- read.csv(fileList[X], row.names = 1, stringsAsFactors = FALSE)
  geneColor <- ifelse(dC$p.adj < 0.1, 'red', 'black')
  if(sID[X] == 'Morphine'){
    gList <- c('Adcy5', 'Foxp1', 'Atp2b1', 'Ppp3ca', 'Rgs7bp', 'Ppp1r1b')
  }
  if(sID[X] %in% c('Aging', '-mtAging')){
    gList <- c('Ndrg4', 'Celf2', 'Pbx1', 'Gria2', 'Malat1', 'Dclk1', 'Ptpro', 'Camk2b', 'Inpp5j', 'Atp1b1', 'Gad1', 'Meis2', 'Ckb', 'Shisa8', 'Nrxn3', 'Prkca', 'Rps19', 'Itm2b')
  }
  if(sID[X] == 'DermalFibroblasts'){
    gList <- c('EEF1A1', 'SERF2', 'B2M')
  }
  if(sID[X] == 'SCC6'){
    gList <- c('H2AFZ', 'DNMT1')
  }
  if(sID[X] == 'NKX21'){
    gList <- NULL
  }
  geneColor[dC$gene %in% gList] <- 'forestgreen'
  genePoint <- (ifelse(dC$p.adj < 0.05, 8, 16))
  geneID <- dC$gene
  geneID[dC$p.adj > 0.05] <- ''
  
  tQ <- -log10(seq(0,1,length.out = nrow(dC)))
  sQ <- -log10(dC$p.value)
  dF <- data.frame(X= tQ, Y= sQ, geneID = geneID)
  #dF <- data.frame(X= sort(rchisq(nrow(dC), df = 1), decreasing = TRUE), Y= dC$FC, geneID = geneID)
  rownames(dF) <- rownames(dC$gene)
  y <- -log10(dC$p.value)
  y <- quantile(y, c(0.2,0.8))
  x <- rev(-log10(qunif(c(0.2,0.8))))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  
  png(paste0('figures/qq',sID[X],'.png'), width = 2835, height = 2025, res = 300, pointsize = 30)
  plotQQ <- ggplot(dF, aes(X,Y, label = geneID)) + 
    geom_point(color = geneColor, pch = genePoint) + 
    theme_bw() + 
    geom_text_repel(segment.color = 'gray60', segment.alpha = 0.5, max.iter = 1e4, aes(fontface = ifelse(dC$gene %in% gList,2,1)), box.padding = .08) + 
    geom_abline(slope = slope, intercept = int, lty = 2) + labs(y=expression(-log[1*0]*" (Observed P-values)"),x=expression(-log[1*0]*" (Expected P-values)"))
  print(plotQQ)
  dev.off()
})
