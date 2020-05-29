library(Matrix)
library(igraph)
library(scTenifoldNet)
dC <- read.csv('results/sym_10x500SCC6_Itensor_Dalignment.csv', row.names = 1, stringsAsFactors = FALSE)
dC <- dC$gene[dC$p.adj < 0.1]

load('datasets/SCC6/SCC6.RData')

X <- as.matrix(O$tensorNetworks$X)
Y <- as.matrix(O$tensorNetworks$Y)

X <- X[!grepl('^RPL|^RPS|^MT-',rownames(X)),!grepl('^RPL|^RPS|^MT-',rownames(X))]
Y <- Y[!grepl('^RPL|^RPS|^MT-',rownames(Y)),!grepl('^RPL|^RPS|^MT-',rownames(Y))]

Y[abs(Y) < quantile(abs(Y),1-(1e-6))] <- 0
X[abs(X) < quantile(abs(X),1-(1e-6))] <- 0

X <- X[rowSums(X) > 0 & colSums(X) > 0,rowSums(X) > 0 & colSums(X) > 0]
Y <- Y[rowSums(Y) > 0 & colSums(Y) > 0,rowSums(Y) > 0 & colSums(Y) > 0]

X <- graph_from_adjacency_matrix(X, weighted = TRUE, diag = FALSE)
Y <- graph_from_adjacency_matrix(Y, weighted = TRUE, diag = FALSE)

intersect(rownames(X[]),rownames(Y[]))

gID <- 'H2AFZ'
sY <- make_ego_graph(Y, nodes = gID, order = 1)[[1]]
sX <- make_ego_graph(X, nodes = gID, order = 1)[[1]]

mX <- readMM('datasets/SCC6/SCC6_CTL.mtx')
rownames(mX) <- readLines('datasets/SCC6/genesSCC6_CTL.txt')

mY <- readMM('datasets/SCC6/SCC6_CTX.mtx')
rownames(mY) <- readLines('datasets/SCC6/genesSCC6_CTX.txt')

mX <- scQC(mX)
mX <- cpmNormalization(mX)

mY <- scQC(mY)
mY <- cpmNormalization(mY)

library(bnlearn)
library(reshape2)

mY <- as.data.frame(as.matrix(t(mY[names(V(sY)),])))
mX <- as.data.frame(as.matrix(t(mX[names(V(sX)),])))

assignDirectionNetworK <- function(igraphNetwork, countMatrix, bootR= 10){
   set.seed(1)
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
sX <- assignDirectionNetworK(igraphNetwork = sX, countMatrix = mX, bootR = 100)

uNet <- igraph::union(sY,sX)
set.seed(1)
lNet <- igraph::layout.graphopt(uNet, charge = 1e-3)
lNet <- t(t(lNet)/apply(abs(lNet),2,max))
rownames(lNet) <- rownames(uNet[])
gY <- names(V(sY))
gX <- names(V(sX))
plot(uNet, layout = lNet, mark.groups = which(gX %in% gY))

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

png('figures/scc6DiffNetworks.png', width = 6000, height = 3000, res = 300, pointsize = 20, bg = NA)
par(mfrow=c(1,2), mar = c(0,0,0,0))
fColor <- ifelse(gX %in% gY, 'darkgoldenrod', NA)
fColor[(gX %in% dC) & fColor == 'darkgoldenrod'] <- 'forestgreen'
plot(sX, layout = lNet[names(V(sX)),], vertex.shape="fcircle", vertex.frame.color=fColor, vertex.frame.width=12, vertex.label = NA, 
     edge.width = NA, edge.arrow.size = 0, rescale = FALSE, xlim = c(min(lNet[,1])*1.05, max(lNet[,1])*1.05), ylim=c(min(lNet[,2])*1.05, max(lNet[,2])*1.05), 
     mark.groups = which(gX %in% gY), mark.col="#C5E5E7", mark.border=NA)
gO2 <- gX
#gO2[gX %in% dC] <- paste0('[', gX[gX %in% dC], ']')
plot(sX, layout = lNet[names(V(sX)),],edge.arrow.size=1,  edge.size = abs(E(sX)$weight)/max(abs(E(sX)$weight)), 
     edge.color = ifelse(E(sX)$weight > 0, 'red', 'blue'), vertex.label=gO2,
     vertex.label.family= 'Arial', vertex.label.color='black', vertex.label.cex= 1, vertex.frame.color = NA, add = TRUE, 
     vertex.color = ifelse(gX %in% dC, 'green3' ,'darkgoldenrod1'), vertex.label.font = ifelse(gX %in% dC, 2,1), 
     rescale = FALSE, xlim = c(min(lNet[,1])*1.05, max(lNet[,1])*1.05), ylim=c(min(lNet[,2])*1.1, max(lNet[,2])*1.05))

fColor <- ifelse(gY %in% gX, 'darkgoldenrod', NA)
fColor[(gY %in% dC) & fColor == 'darkgoldenrod'] <- 'forestgreen'
plot(sY, layout = lNet[names(V(sY)),], vertex.shape="fcircle", vertex.frame.color=fColor, vertex.frame.width=12, vertex.label = NA, edge.width = NA, edge.arrow.size = 0, rescale = FALSE, xlim = c(min(lNet[,1])*1.05, max(lNet[,1])*1.05), ylim=c(min(lNet[,2])*1.05, max(lNet[,2])*1.05), mark.groups = which(gY %in% gX), mark.col="#C5E5E7", mark.border=NA)
gY2 <- gY
#gY2[gY %in% dC] <- paste0('[', gY[gY %in% dC], ']')
plot(sY, layout = lNet[names(V(sY)),], edge.arrow.size=1, edge.size = abs(E(sY)$weight)/max(abs(E(sY)$weight)), edge.color = ifelse(E(sY)$weight > 0, 'red', 'blue'), vertex.label = gY2,
     vertex.label.family= 'Arial', vertex.label.color='black', vertex.label.cex= 1, vertex.frame.color = NA, add = TRUE, vertex.color = ifelse(gY %in% dC, 'green3' ,'darkgoldenrod1'), vertex.label.font = ifelse(gY %in% dC, 2,1), rescale = FALSE, xlim = c(min(lNet[,1])*1.05, max(lNet[,1])*1.05), ylim=c(min(lNet[,2])*1.05, max(lNet[,2])*1.05))


dev.off()

