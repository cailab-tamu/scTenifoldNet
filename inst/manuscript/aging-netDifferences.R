library(Matrix)
library(igraph)
library(scTenifoldNet)
dC <- read.csv('results/nonmit_10X500ANEURON_Itensor_Dalignment.csv', row.names = 1)[,1:30]
dC <- dCoexpression(dC)
dC <- dC$gene[dC$p.adj < 0.1]

Y <- readMM('results/tensorOutput/X_10X500ANEURON_Itensor.mtx')
O <- readMM('results/tensorOutput/Y_10X500ANEURON_Itensor.mtx')
gList <- readLines('results/tensorOutput/genes_10X500ANEURON_Itensor.mtx')
colnames(Y) <- rownames(Y) <- gList
rownames(O) <- colnames(O) <- gList

Y <- as.matrix(Y)
O <- as.matrix(O)

Y[abs(Y) < quantile(abs(Y),1-(1e-6))] <- 0
O[abs(O) < quantile(abs(O),1-(1e-6))] <- 0


Y <- graph_from_adjacency_matrix(Y, weighted = TRUE, diag = FALSE)
O <- graph_from_adjacency_matrix(O, weighted = TRUE, diag = FALSE)

gID <- 'Meis2'
sY <- make_ego_graph(Y, nodes = gID, order = 10)[[1]]
sO <- make_ego_graph(O, nodes = gID, order = 10)[[1]]

mO <- readMM('datasets/agingNeurons/O/oMatrix.mtx')
rownames(mO) <- readLines('datasets/agingNeurons/O/oGenes.txt')


mY <- readMM('datasets/agingNeurons/Y/yMatrix.mtx')
rownames(mY) <- readLines('datasets/agingNeurons/Y/yGenes.txt')

mO <- scQC(mO)
mO <- cpmNormalization(mO)

mY <- scQC(mY)
mY <- cpmNormalization(mY)

library(bnlearn)
library(reshape2)

mY <- as.data.frame(as.matrix(t(mY[names(V(sY)),])))
mO <- as.data.frame(as.matrix(t(mO[names(V(sO)),])))

assignDirectionNetworK <- function(igraphNetwork, countMatrix, bootR= 10){
   bnModel <- boot.strength(countMatrix, algorithm = 'hc', R = bootR)
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
set.seed(1)
lNet <- igraph::layout.graphopt(uNet, charge = 1e-3)
lNet <- t(t(lNet)/apply(abs(lNet),2,max))
rownames(lNet) <- names(V(uNet))
lNet['Uchl1',2] <- lNet['Uchl1',2]*1.1
lNet['Gad1',2] <- lNet['Gad1',2]*1.8
lNet['Itm2b',2] <- lNet['Itm2b',2]*2.2
lNet['Atp1b1',1] <- lNet['Atp1b1',1]*0.2
lNet['Prkca',] <- lNet['Prkca',]*c(0.2,1.8)
lNet['Stmn2',] <- lNet['Stmn2',]*c(1.2,1.6)
lNet['Atp2b1',] <- lNet['Atp2b1',]*c(1.1,1.1)
gY <- names(V(sY))
gO <- names(V(sO))
plot(uNet, layout = lNet, mark.groups = which(gO %in% gY))

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

png('figures/agingDiffNetworks.png', width = 6000, height = 3000, res = 300)
par(mfrow=c(1,2))
plot(sY, layout = lNet[names(V(sY)),], vertex.shape="fcircle", vertex.frame.color=ifelse(gY %in% gO, 'darkgoldenrod', NA), vertex.frame.width=12, vertex.label = NA, edge.width = NA, edge.arrow.size = 0, rescale = FALSE, xlim = c(min(lNet[,1])*1.05, max(lNet[,1])*1.05), ylim=c(min(lNet[,2])*1.05, max(lNet[,2])*1.05), mark.groups = which(gY %in% gO), mark.col="#C5E5E7", mark.border=NA)
gY2 <- gY
gY2[gY %in% dC] <- paste0('[', gY[gY %in% dC], ']')
plot(sY, layout = lNet[names(V(sY)),], edge.arrow.size=1, edge.size = abs(E(sY)$weight)/max(abs(E(sY)$weight)), edge.color = ifelse(E(sY)$weight > 0, 'red', 'blue'), vertex.label = gY2,
     vertex.label.family= 'Arial', vertex.label.color='black', vertex.label.cex= 1, vertex.frame.color = NA, add = TRUE, vertex.color = 'darkgoldenrod1', vertex.label.font = ifelse(gY %in% dC, 2,1), rescale = FALSE, xlim = c(min(lNet[,1])*1.05, max(lNet[,1])*1.05), ylim=c(min(lNet[,2])*1.05, max(lNet[,2])*1.05))
plot(sO, layout = lNet[names(V(sO)),], vertex.shape="fcircle", vertex.frame.color=ifelse(gO %in% gY, 'darkgoldenrod', NA), vertex.frame.width=12, vertex.label = NA, edge.width = NA, edge.arrow.size = 0, rescale = FALSE, xlim = c(min(lNet[,1])*1.05, max(lNet[,1])*1.05), ylim=c(min(lNet[,2])*1.05, max(lNet[,2])*1.05), mark.groups = which(gO %in% gY), mark.col="#C5E5E7", mark.border=NA)
gO2 <- gO
gO2[gO %in% dC] <- paste0('[', gO[gO %in% dC], ']')
plot(sO, layout = lNet[names(V(sO)),],edge.arrow.size=1,  edge.size = abs(E(sO)$weight)/max(abs(E(sO)$weight)), edge.color = ifelse(E(sO)$weight > 0, 'red', 'blue'), vertex.label=gO2,
     vertex.label.family= 'Arial', vertex.label.color='black', vertex.label.cex= 1, vertex.frame.color = NA, add = TRUE, vertex.color = 'darkgoldenrod1', vertex.label.font = ifelse(gO %in% dC, 2,1), rescale = FALSE, xlim = c(min(lNet[,1])*1.05, max(lNet[,1])*1.05), ylim=c(min(lNet[,2])*1.1, max(lNet[,2])*1.05))
dev.off()

