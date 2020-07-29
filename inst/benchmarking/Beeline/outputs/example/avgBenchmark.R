library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(patchwork)

AUPRC <- list.files(pattern = 'AUPRC')
AUPRC <- lapply(AUPRC, function(X){
  D <- read.csv(X)
  colnames(D)[1] <- 'METHOD'
  colnames(D)[2] <- paste0('AUPRC_', colnames(D)[2])
  return(D)
})

D <- AUPRC[[1]]
for(i in seq_along(AUPRC)[-1]){
  D <- merge(D, AUPRC[[i]], by = 'METHOD')
}
AUPRC <- D[order(rowMeans(apply(D[,2:ncol(D)],2,rank)), decreasing = TRUE),]



AUROC <- list.files(pattern = 'AUROC')
AUROC <- lapply(AUROC, function(X){
  D <- read.csv(X)
  colnames(D)[1] <- 'METHOD'
  colnames(D)[2] <- paste0('AUROC_', colnames(D)[2])
  return(D)
})
D <- AUROC[[1]]
for(i in seq_along(AUROC)[-1]){
  D <- merge(D, AUROC[[i]], by = 'METHOD')
}
AUROC <- D[order(rowMeans(apply(D[,2:ncol(D)],2,rank)), decreasing = TRUE),]



TIME <- list.files(pattern = 'Time')
TIME <- lapply(TIME, function(X){
  D <- read.csv(X)
  colnames(D)[1] <- 'METHOD'
  colnames(D)[2] <- paste0('TIME_', colnames(D)[2])
  D[,2] <- 1/D[,2]
  return(D)
})
D <- TIME[[1]]
for(i in seq_along(TIME)[-1]){
  D <- merge(D, TIME[[i]], by = 'METHOD')
}
TIME <- D[order(rowMeans(apply(D[,2:ncol(D)],2,rank)), decreasing = TRUE),]


D <- merge(AUPRC,AUROC, by = 'METHOD')
D <- merge(D, TIME, by = 'METHOD')
rD <- apply(D[,2:ncol(D)],2,rank)
D <- D[order(rowMeans(rD), decreasing = FALSE),]
write.csv(D, file = 'allBenchmark.csv', row.names = FALSE)
timeCol <- which(grepl('TIME', colnames(D)))
D[,timeCol] <- 1/D[,timeCol]
rownames(D) <- D$METHOD
D <- as.data.frame.array(D)
D$METHOD <- factor(D$METHOD, levels = D$METHOD)

X <- data.frame(METHOD = D$METHOD, AUPRC = D$AUPRC_GSD)
X$mCol <- hcl.colors(12, palette = 'Blue-Red 3')[rank(D$AUPRC_GSD)]
X$R <- rank(1/D$AUPRC_GSD)
A <- ggplot(X, aes(AUPRC, METHOD)) + geom_bar(stat="identity", fill= X$mCol) + theme_bw() + 
  theme(axis.title.y = element_blank(), plot.title = element_text(face = 2)) +
  labs(title = 'GSD') +
  geom_text(aes(label=R), hjust=-0.4, cex = 2)

X <- data.frame(METHOD = D$METHOD, AUROC = D$AUROC_GSD)
X$mCol <- hcl.colors(12, palette = 'Blue-Red 3')[rank(D$AUROC_GSD)]
X$R <- rank(1/D$AUROC_GSD)
B <- ggplot(X, aes(AUROC, METHOD)) + geom_bar(stat="identity", fill= X$mCol) + 
  theme_bw() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())  +
  geom_text(aes(label=R), hjust=-0.4, cex=2)

X <- data.frame(METHOD = D$METHOD, TIME = D$TIME_GSD)
X$mCol <- hcl.colors(12, palette = 'Blue-Red 3')[rank(1/D$TIME_GSD)]
X$R <- rank(D$TIME)
C <- ggplot(X, aes(log1p(TIME), METHOD)) + 
  geom_bar(stat="identity", fill= X$mCol) + 
  theme_bw() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) + 
  xlab(expression(log(Seconds + 1))) +
  geom_text(aes(label=R), hjust=-0.1, cex=2)

avgRank <- (rank(1/D$AUPRC_GSD) + rank(1/D$AUROC_GSD) + rank(D$TIME_GSD))/3
X <- data.frame(METHOD = D$METHOD, AVG = avgRank)
X$mCol <- hcl.colors(12, palette = 'Blue-Red 3')[rank(1/X$AVG)]
X$R <- rank(X$AVG)
D <- ggplot(X, aes(AVG, METHOD)) + 
  geom_bar(stat="identity", fill= X$mCol) + 
  theme_bw() + 
  theme(axis.text.y = element_blank(), axis.title.y = element_blank()) + 
  xlab(expression(Average(Ranks))) +
  geom_text(aes(label=R), hjust=-0.1, cex=2)


png('GSD.png', width = 3200,height = 800, res = 300)
A | B | C | D
dev.off()

