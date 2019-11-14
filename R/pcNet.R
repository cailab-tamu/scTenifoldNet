#' @export pcNet
#' @importFrom RSpectra svds
#' @importFrom Matrix Matrix
#' @importFrom pbapply pbsapply
#' @importFrom stats quantile
#' @title pcNet
#' @description Generate a gene regulatory network based on principal components regression.
#' @param X A filtered and normlized gene expression matrix with cells as columns and genes as rows.
#' @param nComp An integer value. The number of principal components in PCA to generate the networks. Should be greater than 2 and lower than the total number of genes.
#' @param scaleScores A boolean value (TRUE/FALSE), if TRUE, the weights will be normalized such that the maximum absolute value is 1.
#' @param symmetric A boolean value (TRUE/FALSE), if TRUE, the weights matrix returned will be symmetric.
#' @param q A decimal value between 0 and 1. Represent the cut-off threshold of top q\% relationships to be returned.
#' @return A gene regulatory network in dgCMatrix format.
#' @references ...
#' @details ...

pcNet <- function(X,
                  nComp = 3,
                  scaleScores = TRUE,
                  symmetric = FALSE,
                  q = 0) {
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
    y <- X[, K]
    Xi <- X
    Xi <- Xi[, -K]
    coeff <- RSpectra::svds(Xi, nComp)$v
    score <- Xi %*% coeff
    score <-
      t(t(score) / (apply(score, 2, function(X) {
        sqrt(sum(X ^ 2))
      }) ^ 2))
    Beta <- colSums(y * score)
    Beta <- coeff %*% (Beta)
    return(Beta)
  }
  X <- (scale(t(X)))
  n <- ncol(X)
  A <- 1 - diag(n)
  B <- pbapply::pbsapply(seq_len(n), pcCoefficients)
  B <- t(B)
  for (K in seq_len(n)) {
    A[K, A[K, ] == 1] = B[K, ]
  }
  if (isTRUE(symmetric)) {
    A <- (A + t(A)) / 2
  }
  absA <- abs(A)
  if (isTRUE(scaleScores)) {
    A <- (A / max(absA))
  }
  A[absA < quantile(absA, q)] <- 0
  diag(A) <- 0
  colnames(A) <- rownames(A) <- gNames
  A <- as(A, 'dgCMatrix')
  return(A)
}
