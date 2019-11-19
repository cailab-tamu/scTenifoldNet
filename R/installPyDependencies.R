#' @export installPyDependencies
#' @importFrom reticulate py_install
#' @title Install \code{Phyton} dependencies
#' @description Install the \code{numpy} and \code{scipy} \code{Phyton} packages using \code{reticulate}.
installPyDependencies <- function(){
  if(!checkPyDependencies()){
    reticulate::py_install(c('numpy', 'scipy'))
  }
}
