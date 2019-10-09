#' @export PCrTdMa
#' @title PCrTdMa
#' @description ...
#' @param X ...
#' @param Y ...
#' @param nNet ...
#' @param nCells ...
#' @param nComp ...
#' @param symmetric ...
#' @param scaleScores ...
#' @param q ...
#' @return ...
#' @references ...
#'
PCrTdMa <- function(X, Y, nNet = 10, nCells = 1000, nComp = 3, symmetric = FALSE, scaleScores = TRUE, q = 0.05){

  xNames <- rownames(X)
  yNames <- rownames(Y)

  sharedGenes <- intersect(xNames, yNames)
  nGenes <- length(sharedGenes)

  X <- X[sharedGenes,]
  Y <- Y[sharedGenes,]

  # X <- as.matrix(X)
  # Y <- as.matrix(Y)

  set.seed(1)
  xList <- makeNetworks(X = X, nCells = nCells, nNet = nNet, nComp = nComp, scaleScores = scaleScores, symmetric = symmetric, q = (1-q))
  yList <- makeNetworks(X = Y, nCells = nCells, nNet = nNet, nComp = nComp, scaleScores = scaleScores, symmetric = symmetric, q = (1-q))

  Tensor <- array(data = 0, dim = c(nGenes, nGenes, 2, nNet))

  for(i in seq_len(nNet)){
    Tensor[,,1,i] <- as.matrix(xList[[i]])
    Tensor[,,2,i] <- as.matrix(yList[[i]])
  }

  Tensor <- rTensor::as.tensor(Tensor)
  Tensor <- rTensor::cp(tnsr = Tensor, num_components = nComp, max_iter = 1000)

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

  tX[abs(tX) < quantile(abs(tX), (1-q))] <- 0
  tY[abs(tY) < quantile(abs(tY), (1-q))] <- 0

  rownames(tX) <- rownames(tY) <- colnames(tX) <- colnames(tY) <- sharedGenes

  mA <- manifoldAlignment(tX, tY, tX+1, tY+1, d = nComp)

  output <- NULL
  output$X <- mA[seq_len(nGenes),]
  rownames(output$X) <- xNames
  output$Y <- mA[-seq_len(nGenes),]
  rownames(output$Y) <- yNames

  return(output)
}
