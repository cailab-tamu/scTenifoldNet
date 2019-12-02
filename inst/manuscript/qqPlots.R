library(Matrix)
library(apcluster)
library(scTenifoldNet)
library(ggplot2)
library(ggrepel)

fileList <- list.files('results/', pattern = '10X500', full.names = TRUE)
sID <- c('Aging', 'DermalFibroblasts', 'Morphine')

sapply(seq_along(fileList), function(X){
  mA <- read.csv(fileList[X], row.names = 1)[,1:30]
  rownames(mA) <- make.unique(toupper(rownames(mA)))
  dC <- dCoexpression(mA, minDist = 1e-5)
  geneColor <- densCols(dC$Z)
  geneColor[dC$p.value < 0.05] <- 'red'
  geneColor <- geneColor
  genePoint <- (ifelse(dC$p.adj < 0.1, 8, 16))
  geneID <- dC$gene
  geneID[dC$p.adj > 0.1] <- ''
  
  tQ <- qqnorm(dC$Z, plot.it = FALSE)
  sQ <- tQ$y
  tQ <- tQ$x
  dF <- data.frame(X= tQ, Y= sQ, geneID = geneID)
  rownames(dF) <- rownames(dC$gene)
  y <- dC$Z
  y <- quantile(y, c(0.25,0.75))
  x <- qnorm(c(0.25,0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  
  png(paste0('figures/qq',sID[X],'.png'), width = 2100, height = 2100, res = 300, pointsize = 12)
  plotQQ <- ggplot(dF, aes(X,Y, label = geneID)) + 
    geom_point(color = geneColor, pch = genePoint) + 
    theme_bw() + 
    geom_text_repel(segment.color = 'gray60', segment.alpha = 0.5, max.iter = 1e4) + 
    geom_abline(slope = slope, intercept = int, lty = 2) + 
    xlab('Theoretical Quantiles (Normal Distribution)') + 
    ylab('Sample Quantiles')
  print(plotQQ)
  dev.off()
  
  
})
