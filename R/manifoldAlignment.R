#' @export manifoldAlignment
#' @importFrom reticulate py_available py_module_available import_from_path
#' @title manifoldAlignment
#' @description ...
#' @param X ...
#' @param Y ...
#' @param wX ...
#' @param wY ...
#' @param d ...
#' @param type ...
#' @return ...
#' @references ...
#' @details ...

manifoldAlignment <- function(X, Y, d=30){
  if(checkPyDependencies()){
    file <- system.file('python/', package = 'scTenifoldNet')
    d <- as.integer(d)
    sharedGenes <- intersect(rownames(X), rownames(Y))
    X <- X[sharedGenes, sharedGenes]
    Y <- Y[sharedGenes, sharedGenes]
    L <- diag(length(sharedGenes))
    # if(type == 'O'){
    #   wX <- X
    #   wY <- Y
    # }
    # if(type == 'D'){
    wX <- X+1
    wY <- Y+1
    # }
    # if(type == 'P'){
    #   wX <- X
    #   wY <- Y
    #   wX[wX != 0] <- wX[wX != 0]+1
    #   wY[wY != 0] <- wY[wY != 0]+1
    # }
    netManifold <- reticulate::import_from_path(module = 'nonLinearManifold', path = file, convert = TRUE)
    alignedNet <- netManifold$nonLinearManifold(X = X, Y = Y, corr = L, num_dims = d, Wx = wX, Wy = wY)
    return(alignedNet)
  }
}
