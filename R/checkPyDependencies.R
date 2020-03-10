checkPyDependencies <- function() {
  if (reticulate::py_available(initialize = TRUE)) {
    pyPackages <- c('numpy', 'scipy', 'functools')
    pyCheck <- sapply(pyPackages, function(X) {
    try(reticulate::py_run_string(paste0('import ', X)))
    })
    pyCheck <- unlist(lapply(pyCheck, function(X){any(class(X) == 'try-error')}))
    if (all(!pyCheck)) {
      return(TRUE)
    } else {
      message(paste0(pyPackages[pyCheck], ' not available'))
    }
  } else {
    stop('Python not available')
  }
}
