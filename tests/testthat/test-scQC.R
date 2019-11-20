test_that("scQC works", {
  # Simulating of a dataset following a negative binomial distribution with high sparcity (~67%)
  nCells = 2000
  nGenes = 100
  set.seed(1)
  X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
  X <- round(X)
  X <- matrix(X, ncol = nCells)
  
  # Warning expected
  expect_warning(scQC(X, minLibSize = 0, removeOutlierCells = FALSE, minPCT = 0, maxMTratio = 1))
  
  # Gene names
  rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
  
  # Input test 1
  expect_equal(class(X), 'matrix')
  X <- scQC(X, minLibSize = 0, removeOutlierCells = FALSE, minPCT = 0, maxMTratio = 1)
  
  # Input test 2
  expect_equal(as.vector(class(X)), 'dgCMatrix')
  X <- scQC(X, minLibSize = 0, removeOutlierCells = FALSE, minPCT = 0, maxMTratio = 1)
  
  # Output test 1 
  X <- scQC(X, minLibSize = 0, removeOutlierCells = FALSE, minPCT = 0, maxMTratio = 1)
  expect_equal(as.vector(class(X)), 'dgCMatrix')
  
  # Output test 3
  expect_equal(dim(X), c(nGenes,nCells))
  
  # Output test 4
  X <- scQC(X, minLibSize = 0, removeOutlierCells = TRUE, minPCT = 0, maxMTratio = 1)
  expect_true(ncol(X) < nCells)
  nCells <- ncol(X)
  
  # Output test 4
  X <- scQC(X, minLibSize = 0, removeOutlierCells = FALSE, minPCT = 0, maxMTratio = 0.1)
  expect_true(ncol(X) < nCells)
  
  # Output test 5
  X <- scQC(X, minLibSize = 0, removeOutlierCells = FALSE, minPCT = 0.5, maxMTratio = 0.1)
  expect_true(ncol(X) < nCells)
  expect_true(nrow(X) < nGenes)
  
})

