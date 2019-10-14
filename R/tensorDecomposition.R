#' @export tensorDecomposition
#' @importFrom rTensor cp
#' 
#' @title tensorDecomposition
#' @description ...
#' @param xList ...
#' @param yList ...
#' @param d ...
#' @return ...
#' @references ...
#' @details ...

tensorDecomposition <- function(xList, yList, d = 5){
  xNets <- length(xList)
  yNets <- length(yList)
  if(xNets != yNets){
    stop('Same number of networks are required in both cases')
  }
  nNet <- unique(c(xNets, yNets))
  xGenes <- unique(unlist(lapply(xList, rownames)))
  yGenes <- unique(unlist(lapply(yList, rownames)))
  sGenes <- intersect(xGenes, yGenes)
  nGenes <- length(sGenes)
  
  Tensor <- array(data = 0, dim = c(nGenes, nGenes, 2, nNet))
  
  for(i in seq_len(nNet)){
    tempX <- matrix(0, nGenes, nGenes)
    rownames(tempX) <- colnames(tempX) <- sGenes
    temp <- as.matrix(xList[[i]])
    tGenes <- sGenes[sGenes %in% rownames(temp)]
    tempX[tGenes,tGenes] <- temp[tGenes,tGenes]
    Tensor[,,1,i] <- tempX
    
    tempY <- matrix(0, nGenes, nGenes)
    rownames(tempY) <- colnames(tempY) <- sGenes
    temp <- as.matrix(yList[[i]])
    tGenes <- sGenes[sGenes %in% rownames(temp)]
    tempY[tGenes,tGenes] <- temp[tGenes,tGenes]
    Tensor[,,2,i] <- tempY
  }
  
  Tensor <- rTensor::as.tensor(Tensor)
  Tensor <- rTensor::cp(tnsr = Tensor, num_components = d, max_iter = 1e3)
  
  tX <- Tensor$est@data[,,1,1]
  tY <- Tensor$est@data[,,2,1]
  
  for(i in seq_len(nNet)[-1]){
    tX <- tX +  Tensor$est@data[,,1,i]
    tY <- tY +  Tensor$est@data[,,2,i]
  }
  
  tX <- tX/nNet
  tY <- tY/nNet
  
  tX <- tX/max(abs(tX))
  tY <- tY/max(abs(tY))
  
  tX <- round(tX,2)
  tY <- round(tY,2)
  
  tX <- as(tX, 'dgCMatrix')
  tY <- as(tY, 'dgCMatrix')
  
  rownames(tX) <- rownames(tY) <- colnames(tX) <- colnames(tY) <- sGenes
  
  tensorOutput <- list()
  tensorOutput$X <- tX
  tensorOutput$Y <- tY
  
  return(tensorOutput)
}

