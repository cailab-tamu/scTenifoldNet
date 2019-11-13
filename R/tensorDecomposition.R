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
  
  # if(type == '4d'){
  #   Tensor <- array(data = 0, dim = c(nGenes, nGenes, 2, nNet))
  #   
  #   for(i in seq_len(nNet)){
  #     tempX <- matrix(0, nGenes, nGenes)
  #     rownames(tempX) <- colnames(tempX) <- sGenes
  #     temp <- as.matrix(xList[[i]])
  #     tGenes <- sGenes[sGenes %in% rownames(temp)]
  #     tempX[tGenes,tGenes] <- temp[tGenes,tGenes]
  #     Tensor[,,1,i] <- tempX
  #     
  #     tempY <- matrix(0, nGenes, nGenes)
  #     rownames(tempY) <- colnames(tempY) <- sGenes
  #     temp <- as.matrix(yList[[i]])
  #     tGenes <- sGenes[sGenes %in% rownames(temp)]
  #     tempY[tGenes,tGenes] <- temp[tGenes,tGenes]
  #     Tensor[,,2,i] <- tempY
  #   }
  #   
  #   Tensor <- rTensor::as.tensor(Tensor)
  #   set.seed(1)
  #   Tensor <- rTensor::cp(tnsr = Tensor, num_components = d, max_iter = 1e3)
  #   
  #   tX <- Tensor$est@data[,,1,1]
  #   tY <- Tensor$est@data[,,2,1]
  #   
  #   for(i in seq_len(nNet)[-1]){
  #     tX <- tX +  Tensor$est@data[,,1,i]
  #     tY <- tY +  Tensor$est@data[,,2,i]
  #   }
  # }
  # 
  # if(type == '3d'){
  #   Tensor <- array(data = 0, dim = c(nGenes, nGenes, 1, 2*nNet))
  #   
  #   for(i in seq_len(nNet)){
  #     tempX <- matrix(0, nGenes, nGenes)
  #     rownames(tempX) <- colnames(tempX) <- sGenes
  #     temp <- as.matrix(xList[[i]])
  #     tGenes <- sGenes[sGenes %in% rownames(temp)]
  #     tempX[tGenes,tGenes] <- temp[tGenes,tGenes]
  #     Tensor[,,1,i] <- tempX
  #     
  #     tempY <- matrix(0, nGenes, nGenes)
  #     rownames(tempY) <- colnames(tempY) <- sGenes
  #     temp <- as.matrix(yList[[i]])
  #     tGenes <- sGenes[sGenes %in% rownames(temp)]
  #     tempY[tGenes,tGenes] <- temp[tGenes,tGenes]
  #     Tensor[,,1,(nNet+i)] <- tempY
  #   }
  #   
  #   Tensor <- rTensor::as.tensor(Tensor)
  #   set.seed(1)
  #   Tensor <- rTensor::cp(tnsr = Tensor, num_components = d, max_iter = 1e3)
  #   
  #   tX <- Tensor$est@data[,,1,1]
  #   tY <- Tensor$est@data[,,1,(nNet+1)]
  #   
  #   for(i in seq_len(nNet)[-1]){
  #     tX <- tX +  Tensor$est@data[,,1,i]
  #     tY <- tY +  Tensor$est@data[,,1,(nNet+i)]
  #   }
  # }
  # 
  # if(type == 'I'){
  tensorX <- array(data = 0, dim = c(nGenes,nGenes,1,nNet))
  tensorY <- array(data = 0, dim = c(nGenes,nGenes,1,nNet))
  
  for(i in seq_len(nNet)){
    tempX <- matrix(0, nGenes, nGenes)
    rownames(tempX) <- colnames(tempX) <- sGenes
    temp <- as.matrix(xList[[i]])
    tGenes <- sGenes[sGenes %in% rownames(temp)]
    tempX[tGenes,tGenes] <- temp[tGenes,tGenes]
    tensorX[,,,i] <- tempX
    
    tempY <- matrix(0, nGenes, nGenes)
    rownames(tempY) <- colnames(tempY) <- sGenes
    temp <- as.matrix(yList[[i]])
    tGenes <- sGenes[sGenes %in% rownames(temp)]
    tempY[tGenes,tGenes] <- temp[tGenes,tGenes]
    tensorY[,,,i] <- tempY
  }
  
  tensorX <- rTensor::as.tensor(tensorX)
  tensorY <- rTensor::as.tensor(tensorY)
  set.seed(1)
  tensorX <- rTensor::cp(tnsr = tensorX, num_components = d, max_iter = 1e3)
  set.seed(1)
  tensorY <- rTensor::cp(tnsr = tensorY, num_components = d, max_iter = 1e3)
  
  tX <- tensorX$est@data[,,,1]
  tY <- tensorY$est@data[,,,1]
  
  for(i in seq_len(nNet)[-1]){
    tX <- tX +  tensorX$est@data[,,,i]
    tY <- tY +  tensorY$est@data[,,,i]
  }
  # }
  
  tX <- tX/nNet
  tY <- tY/nNet
  
  tX <- tX/max(abs(tX))
  tY <- tY/max(abs(tY))
  
  tX <- round(tX,1)
  tY <- round(tY,1)
  
  tX <- as(tX, 'dgCMatrix')
  tY <- as(tY, 'dgCMatrix')
  
  rownames(tX) <- rownames(tY) <- colnames(tX) <- colnames(tY) <- sGenes
  
  tensorOutput <- list()
  tensorOutput$X <- tX
  tensorOutput$Y <- tY
  
  return(tensorOutput)
}

