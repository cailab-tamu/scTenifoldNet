test_that("Rcpp pcNetCoreRcpp matches R compute_gene_coefficients", {

  # Simulated data matching existing test setup
  nCells <- 200
  nGenes <- 100
  nComp  <- 3
  set.seed(42)
  X <- matrix(rnbinom(n = nGenes * nCells, size = 20, prob = 0.98),
              nrow = nGenes, ncol = nCells)
  rownames(X) <- paste0("Gene", seq_len(nGenes))

  # Standardize the same way pcNet does
  X_std <- scale(Matrix::t(X))

  # -----------------------------------------------------------

# Pure R reference: loop over every gene, compute coefficients
  # (mirrors compute_gene_coefficients inside pcNet)
  # -----------------------------------------------------------
  n <- ncol(X_std)
  coeff_R <- matrix(0, nrow = n, ncol = n - 1)

  for (K in seq_len(n)) {
    target_gene   <- X_std[, K]
    design_matrix <- X_std[, -K]

    svd_result         <- RSpectra::svds(design_matrix, nComp)
    principal_components <- svd_result$v

    pc_scores      <- design_matrix %*% principal_components
    score_sq_norms <- colSums(pc_scores^2)
    pc_scores_norm <- pc_scores / rep(score_sq_norms, each = nrow(pc_scores))

    pc_coefficients   <- colSums(target_gene * pc_scores_norm)
    gene_coefficients <- as.vector(principal_components %*% pc_coefficients)

    coeff_R[K, ] <- gene_coefficients
  }

  # -----------------------------------------------------------
  # C++ backend
  # -----------------------------------------------------------
  coeff_Rcpp <- scTenifoldNet:::pcNetCoreRcpp(as.matrix(X_std), nComp)

  # Both must have the same dimensions
  expect_equal(dim(coeff_Rcpp), dim(coeff_R))

  # Values should match within floating-point tolerance
  expect_equal(coeff_Rcpp, coeff_R, tolerance = 1e-6)
})

test_that("pcNet with useRcpp TRUE and FALSE give same network", {
  nCells <- 200
  nGenes <- 100
  set.seed(42)
  X <- matrix(rnbinom(n = nGenes * nCells, size = 20, prob = 0.98),
              nrow = nGenes, ncol = nCells)
  rownames(X) <- paste0("Gene", seq_len(nGenes))

  net_rcpp <- pcNet(X, nComp = 3, scaleScores = FALSE, symmetric = FALSE,
                    q = 0, verbose = FALSE, useRcpp = TRUE)
  net_r    <- pcNet(X, nComp = 3, scaleScores = FALSE, symmetric = FALSE,
                    q = 0, verbose = FALSE, useRcpp = FALSE)

  expect_equal(dim(net_rcpp), dim(net_r))
  expect_equal(as.matrix(net_rcpp), as.matrix(net_r), tolerance = 1e-6)
})
