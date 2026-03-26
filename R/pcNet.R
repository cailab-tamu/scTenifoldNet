#' Principal Component Network Analysis
#'
#' Computes a gene co-expression network using Principal Component Regression.
#' Each gene is regressed against principal components derived from all other genes.
#'
#' @param X A matrix with genes as rows and cells/samples as columns.
#'          Can be a regular matrix or dgCMatrix (sparse).
#'          Must have positive row sums (quality control applied).
#'
#' @param nComp Number of principal components to use for regression.
#'              Must be >= 2 and < number of genes. Default: 3.
#'
#' @param scaleScores If TRUE (default), scales output network by maximum
#'                    absolute value to normalize edge weights to [-1, 1].
#'
#' @param symmetric If TRUE, symmetrizes the network matrix as (A + t(A))/2.
#'                  Default: FALSE.
#'
#' @param q Quantile threshold for filtering edges. Values below this
#'          quantile are set to zero. Range: [0, 1]. Default: 0 (no filtering).
#'
#' @param verbose If TRUE, prints progress updates. Default: FALSE.
#'
#' @param nCores Number of cores for parallelization. If > 1, uses future
#'               backend (multisession). Default: 1 (sequential).
#'
#' @param useRcpp If TRUE (default), uses compiled Rcpp backend for faster
#'                computation. Falls back to pure R if Rcpp backend unavailable.
#'                Default: TRUE.
#'
#' @return A sparse matrix (Matrix::dgCMatrix) representing the gene network.
#'         Entry [i,j] contains the regression coefficient from gene i to gene j.
#'         Diagonal is always zero (no self-loops).
#'
#' @details
#' Algorithm:
#' 1. Standardize input matrix (center and scale by columns)
#' 2. For each gene K:
#'    a. Extract all other genes as design matrix (Xi)
#'    b. Compute truncated SVD of Xi to get nComp principal components
#'    c. Project gene K onto these components
#'    d. Perform OLS regression to get coefficients
#' 3. Assemble coefficients into sparse network matrix
#' 4. Apply optional scaling and filtering
#'
#' @examples
#' \dontrun{
#'   # Create sample data
#'   X <- matrix(rpois(1000*50, 5), nrow=1000, ncol=50)
#'   rownames(X) <- paste0('Gene', 1:1000)
#'
#'   # Compute network
#'   net <- pcNet(X, nComp=3, nCores=4)
#' }
#'
#' @importFrom Matrix rowSums
#' @importFrom Matrix Matrix
#' @importFrom Matrix sparseMatrix
#' @importFrom RSpectra svds
#' @importFrom stats quantile
#' @importFrom furrr future_map
#' @importFrom furrr furrr_options
#' @importFrom cli cli_h1
#' @importFrom cli cli_inform
#' @importFrom cli cli_alert_info
#' @importFrom cli cli_alert_success
#'
#' @export
pcNet <- function(X,
                  nComp = 3,
                  scaleScores = TRUE,
                  symmetric = FALSE,
                  q = 0, 
                  verbose = FALSE,
                  nCores = 1,
                  useRcpp = TRUE) {
  
  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================
  
  # Check quality control: all genes must have at least one count
  if (!all(Matrix::rowSums(X) > 0)) {
    stop("Input matrix contains genes with zero row sums. ",
         "Please apply quality control to remove low-abundance genes.")
  }
  
  # Check input type
  input_class <- class(X)[[1]]
  valid_classes <- c("matrix", "dgCMatrix")
  if (!input_class %in% valid_classes) {
    stop("Input X must be a matrix or dgCMatrix. ",
         "Got: ", input_class)
  }
  
  # Check nComp parameter
  n_genes <- nrow(X)
  if (nComp < 2) {
    stop("nComp must be >= 2. Got: ", nComp)
  }
  if (nComp >= n_genes) {
    stop("nComp must be < number of genes (", n_genes, "). Got: ", nComp)
  }
  
  # ============================================================================
  # DATA PREPARATION
  # ============================================================================
  
  # Store gene names for later
  gene_names <- rownames(X)
  
  # Standardize: transpose to (samples x genes), then center and scale
  X_std <- scale(Matrix::t(X))
  n_genes <- ncol(X_std)  # number of genes after scaling
  
  if (verbose) {
    cli::cli_h1("PCNet - Principal Component Network Analysis")
    cli::cli_inform("Input: {nrow(X_std)} samples x {n_genes} genes")
    cli::cli_inform("Parameters: nComp={nComp}, nCores={nCores}, q={q}")
  }
  
  # ============================================================================
  # PRINCIPAL COMPONENT REGRESSION FUNCTION
  # ============================================================================
  
  # Inner function: compute regression for one gene
  # Input: K = gene index (1 to n_genes)
  # Output: vector of regression coefficients (length n_genes-1)
  compute_gene_coefficients <- function(K) {
    # Target gene to predict
    target_gene <- X_std[, K]
    
    # Design matrix: all genes except target (n_genes-1 features)
    design_matrix <- X_std[, -K]
    
    # Step 1: Compute truncated SVD to get principal components
    # Returns right singular vectors (loadings) for dimension reduction
    svd_result <- RSpectra::svds(design_matrix, nComp)
    principal_components <- svd_result$v
    
    # Step 2: Project design matrix onto principal components
    # Result: (n_samples x nComp) matrix of PC scores
    pc_scores <- design_matrix %*% principal_components
    
    # Normalize PC scores by their squared norms
    # This prevents overfitting to high-variance PCs
    score_sq_norms <- colSums(pc_scores^2)
    pc_scores_normalized <- pc_scores / rep(score_sq_norms, each = nrow(pc_scores))
    
    # Step 3: Regress target gene on normalized PC scores
    # This gives coefficients in the PC space
    pc_coefficients <- colSums(target_gene * pc_scores_normalized)
    
    # Transform PC coefficients back to original gene space
    # Result: (n_genes-1) regression coefficients
    gene_coefficients <- as.vector(principal_components %*% pc_coefficients)
    
    return(gene_coefficients)
  }
  
  # ============================================================================
  # PARALLEL COMPUTATION
  # ============================================================================
  
  if (verbose) {
    cat("Computing regression coefficients...\n")
  }
  
  # Use optimized backend if requested
  use_optimized_backend <- FALSE
  coefficient_matrix <- NULL
  if (isTRUE(useRcpp)) {
    # Try Rcpp first (if compiled)
    result <- tryCatch({
      pcNetCoreRcpp(as.matrix(X_std), nComp)
    }, error = function(e) NULL)
    
    if (!is.null(result)) {
      use_optimized_backend <- TRUE
      coefficient_matrix <- result
      if (verbose) {
        cli::cli_alert_info("Using compiled Rcpp backend")
      }
    }
  }
  
  # Fall back to parallelized R version
  if (!use_optimized_backend) {
    # Set up parallel backend if using multiple cores
    if (nCores > 1) {
      old_plan <- future::plan(future::multisession, workers = nCores)
      on.exit(future::plan(old_plan), add = TRUE)
      if (verbose) {
        cli::cli_alert_info("Using {nCores} cores for parallelization")
      }
    }
    
    if (verbose) {
      cli::cli_alert_info("Computing regression coefficients...")
    }
    
    # Apply regression computation to all genes in parallel
    coefficient_list <- furrr::future_map(
      seq_len(n_genes),
      compute_gene_coefficients,
      .progress = verbose,
      .options = furrr::furrr_options(seed = TRUE)
    )
    
    # Convert list of vectors to matrix (genes x genes-1)
    coefficient_matrix <- do.call(rbind, coefficient_list)
  }
  
  if (verbose) {
    cli::cli_alert_info("Building sparse network matrix...")
  }
  
  # ============================================================================
  # SPARSE MATRIX ASSEMBLY
  # ============================================================================
  
  # Build network matrix using sparse triplet format
  # This is more efficient than creating a dense matrix then converting
  # Pre-allocate triplet indices
  # (each gene contributes n_genes-1 off-diagonal entries)
  n_triplets <- n_genes * (n_genes - 1)
  row_indices <- integer(n_triplets)
  col_indices <- integer(n_triplets)
  values <- numeric(n_triplets)
  
  # Fill triplet format efficiently with direct indexing
  # (avoids expensive vector concatenation in a loop)
  triplet_position <- 0
  for (gene_idx in seq_len(n_genes)) {
    # Off-diagonal column indices for this row
    off_diagonal_cols <- which(seq_len(n_genes) != gene_idx)
    n_off_diag <- length(off_diagonal_cols)
    
    # Direct assignment to pre-allocated vectors
    index_range <- (triplet_position + 1):(triplet_position + n_off_diag)
    row_indices[index_range] <- gene_idx
    col_indices[index_range] <- off_diagonal_cols
    values[index_range] <- coefficient_matrix[gene_idx, ]
    
    triplet_position <- triplet_position + n_off_diag
  }
  
  # Create sparse matrix from triplet representation
  network <- Matrix::sparseMatrix(
    i = row_indices,
    j = col_indices,
    x = values,
    dims = c(n_genes, n_genes)
  )
  
  if (verbose) {
    cli::cli_alert_info("Applying filters and transformations...")
  }
  
  # ============================================================================
  # POST-PROCESSING
  # ============================================================================
  
  # Symmetrize network if requested
  if (isTRUE(symmetric)) {
    network <- (network + Matrix::t(network)) / 2
  }
  
  # Scaling: normalize edge weights to [-1, 1]
  if (isTRUE(scaleScores)) {
    abs_network <- abs(network)
    max_abs_value <- max(abs_network, na.rm = TRUE)
    
    if (is.finite(max_abs_value) && max_abs_value > 0) {
      network <- network / max_abs_value
    }
  }
  
  # Filtering: remove weak edges below quantile threshold
  if (q > 0 && q < 1) {
    abs_network <- abs(network)
    threshold <- stats::quantile(abs_network, q, na.rm = TRUE)
    network[abs_network < threshold] <- 0
  }
  
  # Force diagonal to zero (no self-loops)
  diag(network) <- 0
  
  # Add gene names to rows and columns
  if (!is.null(gene_names)) {
    dimnames(network) <- list(gene_names, gene_names)
  }
  
  if (verbose) {
    cli::cli_alert_success("Done!")
    cli::cli_inform("Network dimensions: {dim(network)[1]} x {dim(network)[2]}")
    cli::cli_inform("Non-zero edges: {length(network@x)}")
  }
  
  # ============================================================================
  # RETURN RESULT
  # ============================================================================
  
  return(network)
}
