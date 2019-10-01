#' @export manifoldAlignment
#' @importFrom reticulate py_available py_module_available source_python
manifoldAlignment <- function(X, Y, d=3){
  if(reticulate::py_available(initialize = TRUE)){
    if(!reticulate::py_module_available('numpy')){
      #reticulate::py_install('numpy')
    } else {
      file <- system.file('python/', package = 'PCrTdMa')
      nonLinearManifold <- reticulate::import_from_path('nonLinearManifold', path = file)
      d <- as.integer(d)
      sharedGenes <- intersect(rownames(X), rownames(Y))
      X <- X[sharedGenes, sharedGenes]
      Y <- Y[sharedGenes, sharedGenes]
      L <- diag(length(sharedGenes))
      nonLinearManifold(X=X, Y=Y, corr=L, num_dims=d, Wx = X, Wy = Y)
    }
  }
}
