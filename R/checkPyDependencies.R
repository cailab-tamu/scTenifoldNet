checkPyDependencies <- function() {
  if (reticulate::py_available(initialize = TRUE)) {
    pyPackages <- c('numpy', 'scipy', 'functools')
    pyCheck <- sapply(pyPackages, function(X) {
      reticulate::py_module_available(X)
    })
    if (all(pyCheck)) {
      return(TRUE)
    } else {
      message(paste0(pyPackages[!pyCheck], ' not available'))
    }
  } else {
    stop('Python not available')
  }
}
