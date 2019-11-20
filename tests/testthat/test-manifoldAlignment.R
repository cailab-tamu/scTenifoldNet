test_that("manifoldAlignment works", {
  # Simulating of a dataset following a negative binomial distribution with high sparcity (~67%)
  nCells = 200
  nGenes = 100
  nNet =  1
  d = 3
  set.seed(1)
  X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
  X <- round(X)
  X <- matrix(X, ncol = nCells)
  rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
  
  X <- cor(t(X))
  mA <- manifoldAlignment(X, X, d = d)
  expect_true(class(mA) == 'matrix')
})
