as.tensor <- function(x,drop=FALSE){
  stopifnot(is.array(x)||is.vector(x))
  tnsr <- list()
  if (is.vector(x)){
    modes <- c(length(x))
    num_modes <- 1L
    tnsr$modes <- modes
    tnsr$num_modes <- num_modes
    tnsr$data <- x
    #new("Tensor", num_modes, modes, data = x)
  }
  else {
    modes <- dim(x)
    num_modes <- length(modes)
    dim1s <- which(modes==1)
    if (drop && (length(dim1s)>0)){
      modes <- modes[-dim1s]
      num_modes <- num_modes-length(dim1s)
      tnsr$modes <- modes
      tnsr$num_modes <- num_modes
      tnsr$data <- array(x,dim=modes)
      #new("list",num_modes,modes,data=array(x,dim=modes))
    } else {
      tnsr$modes <- modes
      tnsr$num_modes <- num_modes
      tnsr$data <- x
    }
  }
  return(tnsr)
}