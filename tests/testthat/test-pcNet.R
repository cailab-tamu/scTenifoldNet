test_that("pCNet works", {
  # Simulating of a dataset following a negative binomial distribution with high sparcity (~67%)
  nCells = 200
  nGenes = 100
  set.seed(1)
  X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
  X <- round(X)
  X <- matrix(X, ncol = nCells)
  rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
  
  # Fake error (a gene with 0 counts)
  eX <- X
  eX[1,] <- 0
  
  # Missing normalization error
  expect_error(pcNet(eX, verbose = FALSE))
  
  # Wrong input error
  expect_error(pcNet(as.data.frame(X), verbose = FALSE))
  
  # Wrong argument error (nComp)
  expect_error(pcNet(X, verbose = FALSE, scaleScores = FALSE, nComp = 1))
  
  # Basic test
  xNet <- pcNet(X, verbose = FALSE, scaleScores = FALSE, symmetric = FALSE, nComp = 2)
  expect_equal(dim(xNet), c(nGenes, nGenes))
  expect_true(all(rownames(xNet) == rownames(X)))
  expect_true(all(colnames(xNet) == rownames(X)))
  
  # Scaling test
  xNet <- pcNet(X, verbose = FALSE, scaleScores = TRUE, symmetric = FALSE, nComp = 2)
  expect_equal(dim(xNet), c(nGenes, nGenes))
  expect_true(all(rownames(xNet) == rownames(X)))
  expect_true(all(colnames(xNet) == rownames(X)))
  expect_true(max(abs(xNet)) == 1)
  
  # Symmetric test
  xNet <- pcNet(X, verbose = FALSE, scaleScores = FALSE, symmetric = TRUE,  nComp = 2)
  expect_equal(dim(xNet), c(nGenes, nGenes))
  expect_true(all(rownames(xNet) == rownames(X)))
  expect_true(all(colnames(xNet) == rownames(X)))
  expect_true(all(xNet[lower.tri(xNet)] == xNet[rev(upper.tri(xNet))]))
  
  # Scaling + Symmetric test
  xNet <- pcNet(X, verbose = FALSE, scaleScores = TRUE, symmetric = TRUE,  nComp = 2)
  expect_equal(dim(xNet), c(nGenes, nGenes))
  expect_true(all(rownames(xNet) == rownames(X)))
  expect_true(all(colnames(xNet) == rownames(X)))
  expect_true(max(abs(xNet)) == 1)
  expect_true(all(xNet[lower.tri(xNet)] == xNet[rev(upper.tri(xNet))]))
  
  # Verbose + Scaling + Symmetric test
  xNet <- pcNet(X, verbose = TRUE, scaleScores = TRUE, symmetric = TRUE,  nComp = 2)
  expect_equal(dim(xNet), c(nGenes, nGenes))
  expect_true(all(rownames(xNet) == rownames(X)))
  expect_true(all(colnames(xNet) == rownames(X)))
  expect_true(max(abs(xNet)) == 1)
  expect_true(all(xNet[lower.tri(xNet)] == xNet[rev(upper.tri(xNet))]))
  
})
