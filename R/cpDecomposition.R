#'@title Canonical Polyadic Decomposition
#'
#'@description Canonical Polyadic (CP) decomposition of a tensor, aka CANDECOMP/PARAFRAC. Approximate a K-Tensor using a sum of \code{num_components} rank-1 K-Tensors. A rank-1 K-Tensor can be written as an outer product of K vectors. There are a total of \code{num_compoents *tnsr@@num_modes} vectors in the output, stored in \code{tnsr@@num_modes} matrices, each with \code{num_components} columns. This is an iterative algorithm, with two possible stopping conditions: either relative error in Frobenius norm has gotten below \code{tol}, or the \code{max_iter} number of iterations has been reached. For more details on CP decomposition, consult Kolda and Bader (2009).
#'@details Uses the Alternating Least Squares (ALS) estimation procedure. A progress bar is included to help monitor operations on large tensors.
#'@name cpDecomposition
#'@importFrom methods is
#'@importFrom stats rnorm
#'@importFrom utils setTxtProgressBar tail txtProgressBar
#'@param tnsr Tensor with K modes
#'@param num_components the number of rank-1 K-Tensors to use in approximation
#'@param max_iter maximum number of iterations if error stays above \code{tol} 
#'@param tol relative Frobenius norm error tolerance
#'@return a list containing the following \describe{
#'\item{\code{lambdas}}{a vector of normalizing constants, one for each component}
#'\item{\code{U}}{a list of matrices - one for each mode - each matrix with \code{num_components} columns}
#'\item{\code{conv}}{whether or not \code{resid} < \code{tol} by the last iteration}
#'\item{\code{norm_percent}}{the percent of Frobenius norm explained by the approximation}
#'\item{\code{est}}{estimate of \code{tnsr} after compression}
#'\item{\code{fnorm_resid}}{the Frobenius norm of the error \code{fnorm(est-tnsr)}}
#'\item{\code{all_resids}}{vector containing the Frobenius norm of error for all the iterations}
#'}
#'@references T. Kolda, B. Bader, "Tensor decomposition and applications". SIAM Applied Mathematics and Applications 2009.

cpDecomposition <- function(tnsr, num_components=NULL,max_iter=25, tol=1e-5){
  kronecker_list <- function(L){
    isvecORmat <- function(x){is.matrix(x) || is.vector(x)}
    stopifnot(all(unlist(lapply(L,isvecORmat))))
    retmat <- L[[1]]
    for(i in 2:length(L)){
      retmat <- kronecker(retmat,L[[i]])
    }
    retmat
  }

  
  fnorm <- function(tnsr){
    arr<-tnsr$data
    sqrt(sum(arr*arr))
  }
  
  rs_unfold <- function(tnsr,m=NULL){
    if(is.null(m)) stop("mode m must be specified")
    num_modes <- tnsr$num_modes
    rs <- m
    cs <- (1:num_modes)[-m]
    unfold(tnsr,row_idx=rs,col_idx=cs)
  }
  
  unfold <- function(tnsr,row_idx=NULL,col_idx=NULL){
    #checks
    rs <- row_idx
    cs <- col_idx
    if(is.null(rs)||is.null(cs)) stop("row and column indices must be specified")
    num_modes <- tnsr$num_modes
    if (length(rs) + length(cs) != num_modes) stop("incorrect number of indices")
    if(any(rs<1) || any(rs>num_modes) || any(cs < 1) || any(cs>num_modes)) stop("illegal indices specified")
    perm <- c(rs,cs)
    if (any(sort(perm,decreasing=TRUE) != num_modes:1)) stop("missing and/or repeated indices")
    modes <- tnsr$modes
    mat <- tnsr$data
    new_modes <- c(prod(modes[rs]),prod(modes[cs]))
    #rearranges into a matrix
    mat <- aperm(mat,perm)
    dim(mat) <- new_modes
    as.tensor(mat)
  }
  
  hadamard_list <- function(L){
    isvecORmat <- function(x){is.matrix(x) || is.vector(x)}
    stopifnot(all(unlist(lapply(L,isvecORmat))))
    retmat <- L[[1]]
    for (i in 2:length(L)){
      retmat <- retmat*L[[i]]
    }
    retmat
  }
  
  khatri_rao_list <- function(L,reverse=FALSE){
    stopifnot(all(unlist(lapply(L,is.matrix))))
    ncols <- unlist(lapply(L,ncol))
    stopifnot(length(unique(ncols))==1)
    ncols <- ncols[1]
    nrows <- unlist(lapply(L,nrow))
    retmat <- matrix(0,nrow=prod(nrows),ncol=ncols)
    if (reverse) L <- rev(L)
    for(j in 1:ncols){
      Lj <- lapply(L,function(x) x[,j])
      retmat[,j] <- kronecker_list(Lj)
    }
    retmat
  }
  
  superdiagonal_tensor <- function(num_modes,len,elements=1L){
    modes <- rep(len,num_modes)
    arr <- array(0, dim = modes)
    if(length(elements)==1) elements <- rep(elements,len)
    for (i in 1:len){
      txt <- paste("arr[",paste(rep("i", num_modes),collapse=","),"] <- ", elements[i],sep="")
      eval(parse(text=txt))
    }
    as.tensor(arr)
  }
  
  ttl<-function(tnsr,list_mat,ms=NULL){
    if(is.null(ms)||!is.vector(ms)) stop ("m modes must be specified as a vector")
    if(length(ms)!=length(list_mat)) stop("m modes length does not match list_mat length")
    num_mats <- length(list_mat)
    if(length(unique(ms))!=num_mats) warning("consider pre-multiplying matrices for the same m for speed")
    mat_nrows <- vector("list", num_mats)
    mat_ncols <- vector("list", num_mats)
    for(i in 1:num_mats){
      mat <- list_mat[[i]]
      m <- ms[i]
      mat_dims <- dim(mat)
      modes_in <- tnsr$modes
      stopifnot(modes_in[m]==mat_dims[2])
      modes_out <- modes_in
      modes_out[m] <- mat_dims[1]
      tnsr_m <- rs_unfold(tnsr,m=m)$data
      retarr_m <- mat%*%tnsr_m
      tnsr <- rs_fold(retarr_m,m=m,modes=modes_out)
    }	
    tnsr	
  }
  
  rs_fold <- function(mat,m=NULL,modes=NULL){
    if(is.null(m)) stop("mode m must be specified")
    if(is.null(modes)) stop("Tensor modes must be specified")
    num_modes <- length(modes)
    rs <- m
    cs <- (1:num_modes)[-m]
    fold(mat,row_idx=rs,col_idx=cs,modes=modes)
  }
  
  fold <- function(mat, row_idx = NULL, col_idx = NULL, modes=NULL){
    #checks
    rs <- row_idx
    cs <- col_idx
    if(is.null(rs)||is.null(cs)) stop("row space and col space indices must be specified")
    if(is.null(modes)) stop("Tensor modes must be specified")
    if(!is(mat,"list")){
      if(!is.matrix(mat))  stop("mat must be of class 'matrix'")
    }else{
      stopifnot(mat$num_modes==2)
      mat <- mat$data			
    }
    num_modes <- length(modes)
    stopifnot(num_modes==length(rs)+length(cs))
    mat_modes <- dim(mat)
    if((mat_modes[1]!=prod(modes[rs])) || (mat_modes[2]!=prod(modes[cs]))) stop("matrix nrow/ncol does not match Tensor modes")
    #rearranges into array
    iperm <- match(1:num_modes,c(rs,cs))
    as.tensor(aperm(array(mat,dim=c(modes[rs],modes[cs])),iperm))
  }
  
  
  if(is.null(num_components)) stop("num_components must be specified")
  stopifnot(is(tnsr,"list"))
  #if (.is_zero_tensor(tnsr)) stop("Zero tensor detected")
  
  #initialization via truncated hosvd
  num_modes <- tnsr$num_modes
  modes <- tnsr$modes
  U_list <- vector("list",num_modes)
  unfolded_mat <- vector("list",num_modes)
  tnsr_norm <- fnorm(tnsr)
  for(m in 1:num_modes){
    unfolded_mat[[m]] <- rs_unfold(tnsr,m=m)$data
    U_list[[m]] <- matrix(rnorm(modes[m]*num_components), nrow=modes[m], ncol=num_components)
  }
  est <- tnsr
  curr_iter <- 1
  converged <- FALSE
  #set up convergence check
  fnorm_resid <- rep(0, max_iter)
  CHECK_CONV <- function(est){
    curr_resid <- fnorm(as.tensor(est$data - tnsr$data))
    fnorm_resid[curr_iter] <<- curr_resid
    if (curr_iter==1) return(FALSE)
    if (abs(curr_resid-fnorm_resid[curr_iter-1])/tnsr_norm < tol) return(TRUE)
    else{ return(FALSE)}
  }	
  #progress bar
  pb <- txtProgressBar(min=0,max=max_iter,style=3)
  #main loop (until convergence or max_iter)
  norm_vec <- function(vec){
    norm(as.matrix(vec))
  }
  while((curr_iter < max_iter) && (!converged)){
    setTxtProgressBar(pb,curr_iter)
    for(m in 1:num_modes){
      V <- hadamard_list(lapply(U_list[-m],function(x) {t(x)%*%x}))
      V_inv <- solve(V)			
      tmp <- unfolded_mat[[m]]%*%khatri_rao_list(U_list[-m],reverse=TRUE)%*%V_inv
      lambdas <- apply(tmp,2,norm_vec)
      U_list[[m]] <- sweep(tmp,2,lambdas,"/")	
      Z <- superdiagonal_tensor(num_modes=num_modes,len=num_components,elements=lambdas)
      est <- ttl(Z,U_list,ms=1:num_modes)
    }
    #checks convergence
    if(CHECK_CONV(est)){
      converged <- TRUE
      setTxtProgressBar(pb,max_iter)
    }else{
      curr_iter <- curr_iter + 1
    }
  }
  if(!converged){setTxtProgressBar(pb,max_iter)}
  close(pb)
  #end of main loop
  #put together return list, and returns
  fnorm_resid <- fnorm_resid[fnorm_resid!=0]
  norm_percent<- (1-(tail(fnorm_resid,1)/tnsr_norm))*100
  invisible(list(lambdas=lambdas, U=U_list, conv=converged, est=est, norm_percent=norm_percent, fnorm_resid = tail(fnorm_resid,1),all_resids=fnorm_resid))
}