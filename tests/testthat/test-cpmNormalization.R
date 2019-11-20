test_that("cpmNormalization works", {
  # Simulating of a dataset following a negative binomial distribution with high sparcity (~67%)
  nCells = 2000
  nGenes = 100
  set.seed(1)
  X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
  X <- round(X)
  X <- matrix(X, ncol = nCells)
  rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
  
  # Input test 1
  expect_equal(class(X), 'matrix')
  X1 <- cpmNormalization(X)
  
  # Output test 1
  expect_equal(as.vector(class(X1)), 'dgCMatrix')
  
  # Output test 2
  expect_equal(dim(X1), c(nGenes, nCells))
  
  # Output test 3
  expect_true(all(round(Matrix::colSums(X1)) == 1e6))
  
  # Input test 2
  X <- as(X, 'dgCMatrix')
  expect_equal(as.vector(class(X)), 'dgCMatrix')
  X2 <- cpmNormalization(X)
  
  # Output test 4
  expect_equal(as.vector(class(X2)), 'dgCMatrix')
  
  # Output test 5
  expect_equal(dim(X2), c(nGenes, nCells))
  
  # Output test 6
  expect_true(all(round(Matrix::colSums(X2)) == 1e6))
  
})
