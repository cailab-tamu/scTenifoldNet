test_that("tensorDecomposition works", {
  # Simulating of a dataset following a negative binomial distribution with high sparcity (~67%)
  nCells = 200
  nGenes = 100
  nNet =  2
  set.seed(1)
  X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
  X <- round(X)
  X <- matrix(X, ncol = nCells)
  rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
  
  X <- makeNetworks(X, nNet = 2, nComp = 2, nCells = 50)
  
  # Wrong number of networks used as input
  expect_error(tensorDecomposition(X, X[1]))
  
  # Basic testing 1 set of networks
  tdOutput <- tensorDecomposition(X)
  expect_true(class(tdOutput) == 'list')
  
  # Basic testing 2 set of networks
  tdOutput <- tensorDecomposition(X, X)
  expect_true(class(tdOutput) == 'list')
})
