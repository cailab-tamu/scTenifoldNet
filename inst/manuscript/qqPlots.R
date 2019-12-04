library(Matrix)
library(apcluster)
library(scTenifoldNet)
library(ggplot2)
library(ggrepel)

fileList <- list.files('results/', pattern = '10X500', full.names = TRUE)
sID <- c('Aging', 'DermalFibroblasts', 'Morphine', '-mtAging')

sapply(seq_along(fileList), function(X){
  mA <- read.csv(fileList[X], row.names = 1)[,1:30]
  rownames(mA) <- make.unique(toupper(rownames(mA)))
  dC <- dCoexpression(mA, minFC = 0.75)
  geneColor <- ifelse(dC$p.value < 0.01, 'red', 'black')
  genePoint <- (ifelse(dC$p.adj < 0.1, 8, 16))
  geneID <- dC$gene
  geneID[dC$p.adj > 0.1] <- ''
  
  tQ <- -log10(seq(0,1,length.out = nrow(dC)))
  sQ <- -log10(dC$p.value)
  dF <- data.frame(X= tQ, Y= sQ, geneID = geneID)
  rownames(dF) <- rownames(dC$gene)
  y <- -log10(dC$p.value)
  y <- quantile(y, c(0.25,0.75))
  x <- rev(-log10(qunif(c(0.25,0.75))))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  
  png(paste0('figures/qq',sID[X],'.png'), width = 2100, height = 1500, res = 300, pointsize = 10)
  plotQQ <- ggplot(dF, aes(X,Y, label = geneID)) + 
    geom_point(color = geneColor, pch = genePoint) + 
    theme_bw() + 
    geom_text_repel(segment.color = 'gray60', segment.alpha = 0.5, max.iter = 1e3) + 
    geom_abline(slope = slope, intercept = int, lty = 2) + labs(y=expression(-log10[1*0]*" (Observed P-values)"),x=expression(-log10[1*0]*" (Expected P-values)"))
  print(plotQQ)
  dev.off()
})
