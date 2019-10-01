#' @export manifoldAlignment
#' @importFrom reticulate py_available py_module_available import_from_path
#'
#' @title manifoldAlignment
#' @description ...
#' @param X ...
#' @param Y ...
#' @param d ...
#' @param method ...
#' @return ...
#' @references ...
#' @details ...

manifoldAlignment <- function(X, Y, d=3, method = 'nonLinear'){
  if(checkPyDependencies()){
    file <- system.file('python/', package = 'PCrTdMa')
    d <- as.integer(d)
    sharedGenes <- intersect(rownames(X), rownames(Y))
    X <- X[sharedGenes, sharedGenes]
    Y <- Y[sharedGenes, sharedGenes]
    L <- diag(length(sharedGenes))
    if(method == 'nonLinear'){
      netManifold <- reticulate::import_from_path(module = 'nonLinearManifold', path = file, convert = TRUE)
      alignedNet <- netManifold$nonLinearManifold(X = X, Y = Y, corr = L, num_dims = d, Wx = X, Wy = Y)
    }
    if(method == 'Linear'){
      netManifold <- reticulate::import_from_path(module = 'linearManifold', path = file, convert = TRUE)
      alignedNet <- netManifold$LinearManifold(X=X, Y=Y, corr = L, num_dims = d, Wx = X, Wy = Y)
    }
    return(alignedNet)
  }
}
