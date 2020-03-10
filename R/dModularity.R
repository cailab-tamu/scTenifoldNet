# library(scTenifoldNet)
# mA <- read.csv('../inst/manuscript/results/10X500DF_Itensor_Dalignment.csv', row.names = 1)[,1:30]
# dR <- dRegulation(mA, minFC = 0)
# 
# dR <- dR[log2(dR$FC) > 1,]
# fMA <- mA[c(paste0('X_', dR$gene), paste0("y_", dR$gene)),1:30]
# X <- fMA
# library(apcluster)
# library(enrichR)
# X <- linKernel(X, normalize = TRUE)
# set.seed(1)
# O <- apcluster(X)#, maxits=1e6, convits = 10)
# C <- O@clusters
# barplot(lengths(C))
# C <- C[lengths(O@clusters) > 10]
# C <- lapply(C, names)
# CX <- lapply(C, function(X){
#   X <- X[grepl('X_',X, ignore.case = TRUE)]
#   gsub('X_','',X)
#   })
# CY <- lapply(C, function(X){
#   X <- X[grepl('Y_',X, ignore.case = TRUE)]
#   gsub('y_','',X, ignore.case = TRUE)
# })
# CY
# 
# JC <- sapply(CX, function(X){
#   sapply(CY, function(Y){
#     A <- table(c(X,Y))
#     mean(A == 2)
#   })
# })
# 
# ComplexHeatmap::Heatmap(JC, row_dend_reorder = FALSE, column_dend_reorder = FALSE)
# 
# JCm <-  reshape2::melt(JC)
# JCm <- JCm[!JCm[,3] %in% c(1,0),]
# 
# L <- unique(c(unique(JCm[,1]), unique(JCm[,2])))
# Ex <- lapply(L, function(i){
#   Ex <- enrichr(CX[[i]], databases = c('BioPlanet_2019', 'KEGG_2019_Human', 'Reactome_2016'))
#   Ex <- do.call(rbind.data.frame, Ex)
#   Ex <- Ex[Ex$Adjusted.P.value < 0.05,]
#   Ex <- Ex[order(Ex$Adjusted.P.value),]
#   Ex$Term
# })
# 
# 
# Ey <- lapply(L, function(i){
#   Ey <- enrichr(CY[[i]], databases = c('BioPlanet_2019', 'KEGG_2019_Human', 'Reactome_2016'))
#   Ey <- do.call(rbind.data.frame, Ey)
#   Ey <- Ey[Ey$Adjusted.P.value < 0.05,]
#   Ey <- Ey[order(Ey$Adjusted.P.value),]
#   Ey$Term
# })
# 
