library(Matrix)
library(igraph)
library(scTenifoldNet)
dC <- read.csv('results/sym10X500morphineNeuron_Itensor_Dalignment.csv', row.names = 1, stringsAsFactors = FALSE)
dC <- dC$gene[dC$p.adj < 0.05]

Y <- readMM('results/tensorOutput/X_10X500morphineNeuron_Itensor.mtx')
O <- readMM('results/tensorOutput/Y_10X500morphineNeuron_Itensor.mtx')
gList <- readLines('results/tensorOutput/genes_10X500morphineNeuron_Itensor.mtx')
colnames(Y) <- rownames(Y) <- gList
rownames(O) <- colnames(O) <- gList

Y <- as.matrix(Y)
O <- as.matrix(O)

Y[abs(Y) < quantile(abs(Y),1-(1e-6))] <- 0
O[abs(O) < quantile(abs(O),1-(1e-6))] <- 0


Y <- graph_from_adjacency_matrix(Y, weighted = TRUE, diag = FALSE)
O <- graph_from_adjacency_matrix(O, weighted = TRUE, diag = FALSE)

Y <- Y[degree(Y) > 0,degree(Y) > 0]
O <- O[degree(O) > 0,degree(O) > 0]

Y <- graph_from_adjacency_matrix(Y, weighted = TRUE, diag = FALSE)
O <- graph_from_adjacency_matrix(O, weighted = TRUE, diag = FALSE)

gID <- 'Ppp3ca'
sY <- make_ego_graph(Y, nodes = gID, order = 1)[[1]]
sO <- make_ego_graph(O, nodes = gID, order = 1)[[1]]

mO <- readMM('datasets/morphineNeurons/morphine/morphine.mtx')
rownames(mO) <- readLines('datasets/morphineNeurons/morphine/morphineGenes.tsv')

mY <- readMM('datasets/morphineNeurons/mock/control.mtx')
rownames(mY) <- readLines('datasets/morphineNeurons/mock/controlGenes.tsv')

mO <- scQC(mO)
mO <- cpmNormalization(mO)

mY <- scQC(mY)
mY <- cpmNormalization(mY)


library(reshape2)

mY <- as.data.frame(as.matrix(t(mY[names(V(sY)),])))
mO <- as.data.frame(as.matrix(t(mO[names(V(sO)),])))

assignDirectionNetworK <- function(igraphNetwork, countMatrix, bootR= 10){
   set.seed(1)
   bnModel <- bnlearn::boot.strength(countMatrix, algorithm = 'hc', R = bootR)
   bnModel$direction <- round(bnModel$direction,1)
   bnModel <- bnModel[bnModel$direction >= 0.5,]
   bY <- graph_from_data_frame(bnModel)
   bY <- as.matrix(bY[])
   igraphNetwork <- as.matrix(igraphNetwork[])
   bY <- bY[rownames(igraphNetwork),colnames(igraphNetwork)]
   igraphNetwork[!(igraphNetwork != 0 & bY != 0)] <- 0
   igraphNetwork <- (graph_from_adjacency_matrix(igraphNetwork, weighted = TRUE))
   return(igraphNetwork)
}

sY <- assignDirectionNetworK(igraphNetwork = sY, countMatrix = mY, bootR = 100)
sO <- assignDirectionNetworK(igraphNetwork = sO, countMatrix = mO, bootR = 100)

uNet <- igraph::union(sY,sO)
set.seed(2)
lNet <- igraph::layout.graphopt(uNet, niter = 1e5)
lNet <- t(t(lNet)/apply(abs(lNet),2,max))
rownames(lNet) <- names(V(uNet))
lNet['Adcy5',] <- lNet['Adcy5',] * c(-1.3,-0.1)
lNet['Rgs7bp',] <- lNet['Rgs7bp',] * c(0.4,0.8)
lNet['Bcl11b',] <- lNet['Bcl11b',] * c(0.5,1.5)
lNet['Vsnl1',] <- lNet['Vsnl1',] * c(1,1.1)
lNet['Ppp1r1b',] <- lNet['Ppp1r1b',] * c(2,1.4)
lNet['Prkcb',] <- lNet['Prkcb',] * c(1,0.7)
lNet['Rasd2',] <- lNet['Rasd2',] * c(0.7,0.8)
plot(uNet, layout = lNet)
gY <- names(V(sY))
gO <- names(V(sO))

library(igraph)

myCircle <- function(coords, v=NULL, params) {
   vertex.color <- params("vertex", "color")
   if (length(vertex.color) != 1 && !is.null(v)) {
      vertex.color <- vertex.color[v]
   }
   vertex.size  <- 1/200 * params("vertex", "size")
   if (length(vertex.size) != 1 && !is.null(v)) {
      vertex.size <- vertex.size[v]
   }
   vertex.frame.color <- params("vertex", "frame.color")
   if (length(vertex.frame.color) != 1 && !is.null(v)) {
      vertex.frame.color <- vertex.frame.color[v]
   }
   vertex.frame.width <- params("vertex", "frame.width")
   if (length(vertex.frame.width) != 1 && !is.null(v)) {
      vertex.frame.width <- vertex.frame.width[v]
   }
   
   mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
          vertex.size, vertex.frame.width,
          FUN=function(x, y, bg, fg, size, lwd) {
             symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                     circles=size, add=TRUE, inches=FALSE)
          })
}

add.vertex.shape("fcircle", clip=igraph.shape.noclip,plot=myCircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))

png('figures/morphineDiffNetworks.png', width = 6000, height = 3000, res = 300, pointsize = 20, bg = NA)
par(mfrow=c(1,2))
fColor <- ifelse(gY %in% gO, 'darkgoldenrod', NA)
fColor[(gY %in% dC) & fColor == 'darkgoldenrod'] <- 'forestgreen'
plot(sY, layout = lNet[names(V(sY)),], vertex.shape="fcircle", vertex.frame.color=fColor, vertex.frame.width=12, vertex.label = NA, edge.width = NA, edge.arrow.size = 0, rescale = FALSE, xlim = c(min(lNet[,1])*1.05, max(lNet[,1])*1.05), ylim=c(min(lNet[,2])*1.05, max(lNet[,2])*1.05), mark.groups = which(gY %in% gO), mark.col="#C5E5E7", mark.border=NA)
gY2 <- gY
#gY2[gY %in% dC] <- paste0('[', gY[gY %in% dC], ']')
plot(sY, layout = lNet[names(V(sY)),], edge.arrow.size=1, edge.size = abs(E(sY)$weight)/max(abs(E(sY)$weight)), edge.color = ifelse(E(sY)$weight > 0, 'red', 'blue'), vertex.label = gY2,
     vertex.label.family= 'Arial', vertex.label.color='black', vertex.label.cex= 1, vertex.frame.color = NA, add = TRUE, vertex.color = ifelse(gY %in% dC, 'green3' ,'darkgoldenrod1'), vertex.label.font = ifelse(gY %in% dC, 2,1), rescale = FALSE, xlim = c(min(lNet[,1])*1.05, max(lNet[,1])*1.05), ylim=c(min(lNet[,2])*1.05, max(lNet[,2])*1.05))
fColor <- ifelse(gO %in% gY, 'darkgoldenrod', NA)
fColor[(gO %in% dC) & fColor == 'darkgoldenrod'] <- 'forestgreen'
plot(sO, layout = lNet[names(V(sO)),], vertex.shape="fcircle", vertex.frame.color=fColor, vertex.frame.width=12, vertex.label = NA, edge.width = NA, edge.arrow.size = 0, rescale = FALSE, xlim = c(min(lNet[,1])*1.05, max(lNet[,1])*1.05), ylim=c(min(lNet[,2])*1.05, max(lNet[,2])*1.05), mark.groups = which(gO %in% gY), mark.col="#C5E5E7", mark.border=NA)
gO2 <- gO
#gO2[gO %in% dC] <- paste0('[', gO[gO %in% dC], ']')
plot(sO, layout = lNet[names(V(sO)),],edge.arrow.size=1,  edge.size = abs(E(sO)$weight)/max(abs(E(sO)$weight)), edge.color = ifelse(E(sO)$weight > 0, 'red', 'blue'), vertex.label=gO2,
     vertex.label.family= 'Arial', vertex.label.color='black', vertex.label.cex= 1, vertex.frame.color = NA, add = TRUE, vertex.color = ifelse(gO %in% dC, 'green3' ,'darkgoldenrod1'), vertex.label.font = ifelse(gO %in% dC, 2,1), rescale = FALSE, xlim = c(min(lNet[,1])*1.05, max(lNet[,1])*1.05), ylim=c(min(lNet[,2])*1.1, max(lNet[,2])*1.05))
dev.off()

