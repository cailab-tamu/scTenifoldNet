test_that("makeNetworks works", {
  # Simulating of a dataset following a negative binomial distribution with high sparcity (~67%)
  nCells = 200
  nGenes = 100
  nNet =  2
  set.seed(1)
  X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
  X <- round(X)
  X <- matrix(X, ncol = nCells)
  
  # Missing gene ids error
  expect_error(makeNetworks(X))
  
  rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
  
  # Wrong argument error (nComp)
  expect_error(makeNetworks(X, nComp = 1))
  
  # Basic tests
  xNet <- makeNetworks(X, nNet = nNet, nCells = 50)
  expect_true(class(xNet) == 'list')
  expect_true(length(xNet) == nNet)
  dimTest <- lapply(xNet, function(C){
    all(dim(C) == c(nGenes, nGenes))
  })
  dimTest <- unlist(dimTest)
  expect_true(all(dimTest))
})
