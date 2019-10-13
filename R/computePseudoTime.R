#' #' @importClassesFrom Biobase AnnotatedDataFrame
#' #' @importFrom Biobase pData
#' #' @importFrom monocle newCellDataSet  reduceDimension orderCells
#' #' @importFrom BiocGenerics estimateSizeFactors
#' #' @importFrom methods new
#' #' @importFrom VGAM negbinomial.size
#' 
#' computePseudoTime <- function(X){
#'   if(is.null(rownames(X))){
#'     rownames(X) <- paste0('gene_',seq_len(nrow(X)))
#'   }
#'   X <- X[apply(X!=0, 1, sum) > 0,]
#'   fd <- data.frame('gene_short_name' = rownames(X))
#'   rownames(fd) <- rownames(X)
#'   fd <- methods::new("AnnotatedDataFrame", data = fd)
#'   cds <- monocle::newCellDataSet(as.matrix(X), featureData = fd, expressionFamily = VGAM::negbinomial.size())
#'   cds <- BiocGenerics::estimateSizeFactors(cds)
#'   cds <- monocle::reduceDimension(cds, reduction_method = "DDRTree", verbose = TRUE, max_components = 2)
#'   cds <- monocle::orderCells(cds)
#'   pseudoT <- Biobase::pData(cds)
#'   pseudoT <- cbind(pseudoT,t(cds@reducedDimS[,rownames(pseudoT)]))
#'   colnames(pseudoT) <- c("sizeFactor", "pseudoTime", 'stateC', 'DDRT1','DDRT2')
#'   return(pseudoT)
#' }
