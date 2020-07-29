library(GGally)
library(ggplot2)
library(ggrepel)
library(patchwork)

inputFile <- list.files(pattern = '.RData')

inputData <- lapply(inputFile, function(X){
  load(X)
  O <- O$diffRegulation
  Z <- O$Z
  names(Z) <- O$gene
  return(Z)
})

gNames <- table(unlist(lapply(inputData, names)))
gNames <- gNames[gNames %in% length(inputData)]
gNames <- names(gNames)

inputData <- lapply(inputData, function(X){X[gNames]})
inputData <- as.data.frame(inputData)
colnames(inputData) <- gsub('.RData','',inputFile)

PC <- prcomp(t(apply(inputData,2,rank)))$x
PC <- as.data.frame(PC[,1:2])
PC$ID <- rownames(PC)

png('K.png', width = 400*8, height = 400*8, res = 300)
ggpairs(inputData, upper = list(continuous = "cor", corMethod = "spearman"))
dev.off()

png('PC.png', width = 1100, height = 1000, res = 200)
ggplot(PC, aes(PC1,PC2,label=ID)) + geom_point() + geom_text_repel() + theme_bw() + labs(title = 'K tensor') + theme(plot.title = element_text(face = 2))
dev.off()
