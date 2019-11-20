test_that("multiplication works", {
  # Simulating of a dataset following a negative binomial distribution with high sparcity (~67%)
  nCells = 200
  nGenes = 100
  set.seed(1)
  X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
  X <- round(X)
  X <- matrix(X, ncol = nCells)
  rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
  
  out <- scTenifoldNet(X, X, qc_minLibSize = 30, nc_nNet = 2, nc_nCells = 50, dc_minDist = 0)
  expect_equal(class(out), 'list')
})
