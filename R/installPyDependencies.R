#' @export installPyDependencies
#' @importFrom reticulate py_install
#' 
installPyDependencies <- function(){
  if(!checkPyDependencies()){
    reticulate::py_install(c('numpy', 'scipy', 'functools')) 
  }
}