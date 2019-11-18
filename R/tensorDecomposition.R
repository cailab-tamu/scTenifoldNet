#' @export tensorDecomposition
#' @importFrom rTensor cp
#'
#' @title CANDECOMP/PARAFAC (CP) Tensor Decomposition of Gene Regulatory Networks.
#' @description Generate weight-averaged denoised gene regulatory networks using CANDECOMP/PARAFAC (CP) Tensor Decomposition.The \code{tensorDecomposition} function takes one or two lists of gene regulatory matrices, if two list are provided, the shared genes are selected and the CP tensor decomposition is performed independently for each list (3d-tensor). The tensor decomposed matrices are then averaged to generate weight-averaged denoised networks.
#' @param xList A list of gene regulatory networks.
#' @param yList Optional. A list of gene regulatory networks.
#' @param K The number of rank-one tensors used to approximate the data using CANDECOMP/PARAFAC (CP) Tensor Decomposition,
#' @param maxError A decimal value between 0 and 1. Defines the relative Frobenius norm error tolerance
#' @param maxIter An integer value. Defines the maximum number of iterations if error stay above \code{maxError}.
#' @return A list of weight-averaged denoised gene regulatory networks.
#' @references 
#' \itemize{
#' \item Kolda, Tamara G., and Brett W. Bader. "Tensor decompositions and applications." SIAM review 51.3 (2009): 455-500.
#' \item Morup, Morten. "Applications of tensor (multiway array) factorizations and decompositions in data mining." Wiley Interdisciplinary Reviews: Data Mining and Knowledge Discovery 1.1 (2011): 24-40.
#' }
#' @details CANDECOMP/PARAFRAC (CP) tensor decomposition approximate a K-Tensor using a sum of \code{d} rank-1 K-Tensors. A rank-1 K-Tensor can be written as an outer product of K vectors. This is an iterative algorithm, with two possible stopping conditions: either relative error in Frobenius norm has gotten below \code{maxError}, or the \code{maxIter} number of iterations has been reached. For more details on CP decomposition, consult Kolda and Bader (2009) and Morup (2011).
#' 
tensorDecomposition <- function(xList, yList = NULL, K = 5, maxError = 1e-5, maxIter = 1e3){
  xNets <- length(xList)
  if(!is.null(yList)){
    yNets <- length(yList)
    if(xNets != yNets){
      stop('Same number of networks are required in both cases')
    }
    nNet <- unique(c(xNets, yNets))
    xGenes <- unique(unlist(lapply(xList, rownames)))
    yGenes <- unique(unlist(lapply(yList, rownames)))
    sGenes <- intersect(xGenes, yGenes)
  } else {
    nNet <- xNets
    xGenes <- unique(unlist(lapply(xList, rownames)))
    sGenes <- xGenes
  }
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
  #   Tensor <- rTensor::cp(tnsr = Tensor, num_components = K, max_iter = maxIter, tol = maxError)
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
  #   Tensor <- rTensor::cp(tnsr = Tensor, num_components = K, max_iter = maxIter, tol = maxError)
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
  if(!is.null(yList)){
    tensorY <- array(data = 0, dim = c(nGenes,nGenes,1,nNet))
  }
  
  for(i in seq_len(nNet)){
    tempX <- matrix(0, nGenes, nGenes)
    rownames(tempX) <- colnames(tempX) <- sGenes
    temp <- as.matrix(xList[[i]])
    tGenes <- sGenes[sGenes %in% rownames(temp)]
    tempX[tGenes,tGenes] <- temp[tGenes,tGenes]
    tensorX[,,,i] <- tempX
    
    if(!is.null(yList)){
      tempY <- matrix(0, nGenes, nGenes)
      rownames(tempY) <- colnames(tempY) <- sGenes
      temp <- as.matrix(yList[[i]])
      tGenes <- sGenes[sGenes %in% rownames(temp)]
      tempY[tGenes,tGenes] <- temp[tGenes,tGenes]
      tensorY[,,,i] <- tempY
    }
    
  }
  
  tensorX <- rTensor::as.tensor(tensorX)
  set.seed(1)
  tensorX <- rTensor::cp(tnsr = tensorX, num_components = K, max_iter = maxIter, tol = maxError)
  tX <- tensorX$est@data[,,,1]
  for(i in seq_len(nNet)[-1]){
    tX <- tX +  tensorX$est@data[,,,i]
  }
  tX <- tX/nNet
  tX <- tX/max(abs(tX))
  tX <- round(tX,1)
  tX <- as(tX, 'dgCMatrix')
  rownames(tX) <- colnames(tX) <- sGenes
  
  if(!is.null(yList)){
    tensorY <- rTensor::as.tensor(tensorY)
    set.seed(1)
    tensorY <- rTensor::cp(tnsr = tensorY, num_components = K, max_iter = 1e3)
    tY <- tensorY$est@data[,,,1]
    for(i in seq_len(nNet)[-1]){
      tY <- tY +  tensorY$est@data[,,,i]
    }
    tY <- tY/nNet
    tY <- tY/max(abs(tY))
    tY <- round(tY,1)
    tY <- as(tY, 'dgCMatrix')  
    rownames(tY) <- colnames(tY) <- sGenes
  }
  
  tensorOutput <- list()
  tensorOutput$X <- tX
  if(!is.null(yList)){
    tensorOutput$Y <- tY
  }
  return(tensorOutput)
}

