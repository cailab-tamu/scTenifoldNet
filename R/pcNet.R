#' @export pcNet
#' @importFrom RSpectra svds
#' @importFrom Matrix Matrix
#' @importFrom pbapply pbsapply
#' @importFrom stats quantile
#' @title Computes a gene regulatory network based on principal component regression
#' @description This function computes a gene regulatory network based on principal component regression (PCR), a technique based on principal component analysis. In PCR, the outcome variable is regressed over a \code{nComp} number of for principal components computed from a set of covariates to estimate the unknown regression coefficients in the model. \code{pcNet} function computes the PCR coefficients for each gene one at a time using all the others as covariates, to construct an all by all gene regulatory network.
#' @param X A filtered and normalized gene expression matrix with cells as columns and genes as rows.
#' @param nComp An integer value. The number of principal components in PCA to generate the networks. Should be greater than 2 and lower than the total number of genes.
#' @param scaleScores A boolean value (\code{TRUE/FALSE}), if TRUE, the weights will be normalized such that the maximum absolute value is 1.
#' @param symmetric A boolean value (TRUE/FALSE), if TRUE, the weights matrix returned will be symmetric.
#' @param q A decimal value between 0 and 1. Defines the cut-off threshold of top q\% relationships to be returned.
#' @param verbose A boolean value (TRUE/FALSE), if TRUE, a progress bar is shown.
#' @return A gene regulatory network in dgCMatrix format.
#' @references 
#' \itemize{
#' \item Jolliffe, Ian T. "A note on the use of principal components in regression." Journal of the Royal Statistical Society: Series C (Applied Statistics) 31.3 (1982): 300-303.
#' }
#' @details Principal component regression may be broadly divided into three major steps: \enumerate{
#' \item Perform PCA on the observed covariates data matrix to obtain \code{nComp} number of the principal components.
#' \item Regress the observed vector of outcomes on the selected principal components as covariates, using ordinary least squares regression to get a vector of estimated regression coefficients
#' \item Transform this vector back to the scale of the actual covariates, using the eigenvectors corresponding to the selected principal components to get the final PCR estimator for estimating the regression coefficients characterizing the original model.
#' }

pcNet <- function(X,
                  nComp = 3,
                  scaleScores = TRUE,
                  symmetric = FALSE,
                  q = 0, verbose = TRUE) {
  if (!all(rowSums(X) > 0)) {
    stop('Quality control has not been applied over the matrix.')
  }
  if (!is.matrix(X)) {
    stop('Input should be a matrix with cells as columns and genes as rows')
  }
  if (nComp < 2 | nComp >= nrow(X)) {
    stop('nCom should be greater than 2 and lower than the total number of genes')
  }
  gNames <- rownames(X)
  pcCoefficients <- function(K) {
    # Taking out the gene to be regressed out
    y <- X[, K]
    Xi <- X
    Xi <- Xi[, -K]
    # Step 1: Perform PCA on the observed covariates data matrix to obtain $n$ number of the principal components.
    coeff <- RSpectra::svds(Xi, nComp)$v
    score <- Xi %*% coeff
    # Step 2: Regress the observed vector of outcomes on the selected principal components as covariates, using ordinary least squares regression to get a vector of estimated regression coefficients.
    score <-
      t(t(score) / (apply(score, 2, function(X) {
        sqrt(sum(X ^ 2))
      }) ^ 2))
    # Step 3: Transform this vector back to the scale of the actual covariates, using the eigenvectors corresponding to the selected principal components to get the final PCR estimator for estimating the regression coefficients characterizing the original model.
    Beta <- colSums(y * score)
    Beta <- coeff %*% (Beta)
    
    return(Beta)
  }
  
  # Standardizing the data
  X <- (scale(t(X)))
  
  # Identify the number of rows in the input matrix
  n <- ncol(X)
  
  # Generate the output matrix
  A <- 1 - diag(n)
  
  # Apply the principal component regression for each gene
  if(verbose){
    B <- pbapply::pbsapply(seq_len(n), pcCoefficients)  
  } else {
    B <- sapply(seq_len(n), pcCoefficients)  
  }
  
  # Transposition of the Beta coefficient matrix
  B <- t(B)
  
  # Replacing the values in the output matrix
  for (K in seq_len(n)) {
    A[K, A[K, ] == 1] = B[K, ]
  }
  
  # Making the output matrix symmetric
  if (isTRUE(symmetric)) {
    A <- (A + t(A)) / 2
  }
  
  # Absolute values for scaling and filtering
  absA <- abs(A)
  
  # Scaling the output matrix
  if (isTRUE(scaleScores)) {
    A <- (A / max(absA))
  }
  
  # Filtering the output matrix
  A[absA < quantile(absA, q)] <- 0
  
  # Setting the diagonal to be 0
  diag(A) <- 0
  
  # Adding names
  colnames(A) <- rownames(A) <- gNames
  
  # Making the output a sparse matrix
  A <- as(A, 'dgCMatrix')
  
  # Return
  return(A)
}
