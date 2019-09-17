#' @export pcNet
#' @importFrom RSpectra svds
#' @importFrom Matrix Matrix
#' @importFrom pbapply pbsapply
#' @title pcNet
#' @description ...
#' @param X ...
#' @param nCom ...
#' @param scaleScores ...
#' @param symmetric ...
#' @param q ...
#' @return ...
#' @references ...
#' @details ...
#' @examples
#' # Generating input matrix
#' set.seed(1)
#' inputMatrix <- matrix(data = rnbinom(n = 1000, size = 10, prob = .9), nrow = 20)
#' rownames(inputMatrix) <- paste0('Gene_', seq_len(nrow(inputMatrix)))
#'
#' # Getting the adjacency matrix
#' pcNet(inputMatrix)
#' # Getting the adjacency matrix filtered to the top 1%
#' pcNet(inputMatrix, q = 0.99)
#' # Getting the adjacency matrix filtered to the top 5%
#' pcNet(inputMatrix, q = 0.95)
#' # Getting the adjacency matrix filtered to the top 10%
#' pcNet(inputMatrix, q = 0.90)

pcNet <- function(X,
                  nCom = 3,
                  scaleScores = TRUE,
                  symmetric = FALSE,
                  q = 0) {
  # Input as matrix
  if (!is.matrix(X)) {
    stop('Input should be a matrix of n x m where n are genes and m are cells')
  }
  if (nCom < 2 | nCom >= nrow(X)) {
    stop('nCom should be greater than 2 and lower than the number of genes')
  }
  gNames <- rownames(X)
  pcCoefficients <- function(K) {
    y <- X[, K]
    Xi <- X
    Xi <- Xi[,-K]
    coeff <- RSpectra::svds(Xi, nCom)$v
    score <- Xi %*% coeff
    score <-
      t(t(score) / (apply(score, 2, function(X) {
        sqrt(sum(X ^ 2))
      }) ^ 2))
    Beta <- colSums(y * score)
    Beta <- coeff %*% (Beta)
    return(Beta)
  }
  n <- ncol(X)
  X <- (scale(t(X)))
  A <- 1 - diag(n)
  B <- pbapply::pbsapply(seq_len(n), pcCoefficients)
  B <- t(B)
  for (K in seq_len(n)) {
    A[K, A[K,] == 1] = B[K,]
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
  A <- Matrix::Matrix(A, sparse = TRUE)
  return(A)
}
