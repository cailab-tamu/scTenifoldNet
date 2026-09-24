test_that("pCNet works", {
  # Simulating of a dataset following a negative binomial distribution with high sparcity (~67%)
  nCells = 200
  nGenes = 100
  set.seed(1)
  X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
  X <- round(X)
  X <- matrix(X, ncol = nCells)
  rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
  
  # A gene with 0 counts gets no edges, and the rest of the network is unchanged
  eX <- X
  eX[1,] <- 0
  eNet <- pcNet(eX, verbose = FALSE, scaleScores = FALSE, nComp = 2)
  expect_equal(dim(eNet), c(nGenes, nGenes))
  expect_true(all(eNet[1, ] == 0) && all(eNet[, 1] == 0))
  expect_equal(as.matrix(eNet[-1, -1]),
               as.matrix(pcNet(X[-1, ], verbose = FALSE, scaleScores = FALSE, nComp = 2)))
  
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

test_that("pcNet matches per-gene SVD principal component regression", {
  nCells <- 200
  nGenes <- 100
  nComp  <- 3
  set.seed(42)
  X <- matrix(rnbinom(n = nGenes * nCells, size = 20, prob = 0.98),
              nrow = nGenes, ncol = nCells)
  rownames(X) <- paste0("Gene", seq_len(nGenes))

  # Reference: regress each gene on the top principal components of the others
  X_std <- scale(Matrix::t(X))
  reference <- matrix(0, nGenes, nGenes)
  for (K in seq_len(nGenes)) {
    design_matrix <- X_std[, -K]
    principal_components <- svd(design_matrix, nu = 0, nv = nComp)$v
    pc_scores <- design_matrix %*% principal_components
    pc_coefficients <- colSums(X_std[, K] * pc_scores) / colSums(pc_scores^2)
    reference[K, -K] <- principal_components %*% pc_coefficients
  }

  xNet <- pcNet(X, nComp = nComp, scaleScores = FALSE, verbose = FALSE)
  expect_equal(unname(as.matrix(xNet)), reference, tolerance = 1e-8)

  # Parallel blocks give the same network
  expect_equal(pcNet(X, nComp = nComp, verbose = FALSE, nCores = 2),
               pcNet(X, nComp = nComp, verbose = FALSE))
})
