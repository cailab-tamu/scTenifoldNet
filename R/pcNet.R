#' Principal Component Network Analysis
#'
#' Builds a gene co-expression network with principal component regression.
#' Every gene is regressed on the leading principal components of all the
#' other genes, and the regression coefficients become the edge weights.
#'
#' @param X Expression matrix with genes in rows and cells (or samples) in
#'   columns. Either a base `matrix` or a sparse `dgCMatrix`. Genes with the
#'   same value in every cell (for example, all zeros) carry no information:
#'   they are left out of the regressions and get no edges.
#' @param nComp Number of principal components used in each regression.
#'   Must be at least 2 and smaller than the number of non-constant genes.
#'   Default: 3.
#' @param scaleScores If `TRUE` (default), divides the network by its largest
#'   absolute value so that edge weights lie in \[-1, 1\].
#' @param symmetric If `TRUE`, replaces the network `A` by `(A + t(A)) / 2`.
#'   Default: `FALSE`.
#' @param q Quantile used to drop weak edges. Edges whose absolute weight is
#'   below this quantile are set to zero. Must be in \[0, 1\]; values of 0 or 1
#'   disable filtering. Default: 0.
#' @param priorNetwork Optional prior gene regulatory network: a data.frame
#'   with two columns, regulators and targets. When `q` filters edges, edges
#'   in the prior network keep their weight. Default: `NULL`.
#' @param verbose If `TRUE`, prints progress messages. Default: `FALSE`.
#' @param nCores Number of cores. With more than one core, the genes are split
#'   into blocks that are processed in parallel (forked processes through
#'   `parallel::mclapply`, or a PSOCK cluster on Windows). If R uses a
#'   multithreaded BLAS (OpenBLAS, Accelerate, MKL), keep `nCores = 1` and let
#'   the BLAS use the cores instead: combining both runs too many threads.
#'   Default: 1.
#' @param useRcpp Ignored. Kept so that existing calls keep working.
#'
#' @return A sparse `dgCMatrix` with all input genes in rows and columns.
#'   Entry `[i, j]` is the coefficient of gene `j` in the regression of gene
#'   `i`. The diagonal is always zero, and so are the rows and columns of
#'   constant genes.
#'
#' @details
#' **What is computed.** The expression matrix is standardized so that each
#' gene has mean 0 and standard deviation 1. For each gene `k`, let `y` be its
#' expression and `A` the matrix of all other genes. The coefficients are
#'
#' `beta_k = V diag(1 / s) U' y`
#'
#' where `U`, `s` and `V` are the leading `nComp` singular vectors and values
#' of `A`. This is ordinary least squares on the top principal components,
#' mapped back to gene space.
#'
#' **How it is computed.** Running one SVD per gene is slow. Instead, note
#' that removing gene `k` changes the cell-by-cell Gram matrix by a single
#' rank-one term: `A A' = X X' - y y'`. The function therefore:
#'
#' 1. Eigendecomposes one Gram matrix, `X X' = Q D Q'`, working in whichever
#'    of the cell or gene dimensions is smaller.
#' 2. For each gene, finds the top `nComp` eigenvalues of `D - z z'`, where
#'    `z = Q' y`. These are the roots of a scalar equation (the "secular
#'    equation"), solved for all genes at once. The matching eigenvectors
#'    have a closed form, so no per-gene SVD is needed.
#' 3. Obtains all coefficients from a single matrix product.
#'
#' The result equals the per-gene SVD approach up to floating point error.
#'
#' **Constant genes.** A gene with the same value in every cell cannot be
#' standardized (its standard deviation is 0). Such genes are removed before
#' the steps above, including scaling and quantile filtering, and are added
#' back at the end as rows and columns of zeros. The network among the other
#' genes is therefore the same as if the constant genes had been removed from
#' `X` beforehand.
#'
#' @examples
#' \dontrun{
#'   X <- matrix(rpois(1000 * 50, 5), nrow = 1000, ncol = 50)
#'   rownames(X) <- paste0("Gene", 1:1000)
#'   net <- pcNet(X, nComp = 3)
#' }
#'
#' @importFrom Matrix rowSums sparseMatrix t
#' @importFrom methods new as
#' @importFrom stats quantile
#' @importFrom parallel mclapply makePSOCKcluster parLapply stopCluster
#' @importFrom cli cli_h1 cli_inform cli_alert_info cli_alert_success cli_abort
#'
#' @export
pcNet <- function(X,
                  nComp = 3,
                  scaleScores = TRUE,
                  symmetric = FALSE,
                  q = 0,
                  priorNetwork = NULL,
                  verbose = FALSE,
                  nCores = 1,
                  useRcpp = TRUE) {

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------

  input_class <- class(X)[[1]]
  if (!input_class %in% c("matrix", "dgCMatrix")) {
    stop("Input X must be a matrix or dgCMatrix. Got: ", input_class)
  }

  n_genes <- nrow(X)
  gene_names <- rownames(X)

  if (nComp < 2) {
    stop("nComp must be >= 2. Got: ", nComp)
  }

  # Constant genes (e.g. all zeros) are set aside; everything up to the final
  # step works on the remaining genes only
  standardized <- .pcNetStandardize(X)
  X_std <- standardized$X_std
  used <- standardized$non_constant
  n_used <- sum(used)
  rm(standardized)

  if (nComp >= n_used) {
    stop("nComp must be < number of non-constant genes (", n_used, "). Got: ", nComp)
  }

  prior_matrix <- NULL
  if (!is.null(priorNetwork)) {
    prior_matrix <- .pcNetPriorMatrix(priorNetwork, gene_names)[used, used, drop = FALSE]
  }

  nCores <- max(1L, min(as.integer(nCores), n_used))

  # ---------------------------------------------------------------------------
  # Regression coefficients
  # ---------------------------------------------------------------------------

  if (verbose) {
    cli::cli_h1("PCNet - Principal Component Network Analysis")
    cli::cli_inform("Input: {nrow(X_std)} samples x {n_genes} genes")
    if (n_used < n_genes) {
      cli::cli_alert_info("{n_genes - n_used} constant gene{?s} (e.g. all zeros) will have no edges")
    }
    cli::cli_inform("Parameters: nComp={nComp}, nCores={nCores}, q={q}")
    cli::cli_alert_info("Eigendecomposing Gram matrix...")
  }

  basis <- .pcNetEigenBasis(X_std)
  rm(X_std)

  if (verbose) {
    cli::cli_alert_info("Computing regression coefficients...")
    if (nCores > 1L) {
      cli::cli_alert_info("Using {nCores} cores for parallelization")
    }
  }

  network <- .pcNetCoefficients(basis$Z, basis$d, nComp, nCores)
  rm(basis)
  diag(network) <- 0

  # ---------------------------------------------------------------------------
  # Post-processing
  # ---------------------------------------------------------------------------
  # Done on the dense matrix. The values match what the equivalent sparse
  # operations would give.

  if (verbose) {
    cli::cli_alert_info("Applying filters and transformations...")
  }

  if (isTRUE(symmetric)) {
    network <- (network + t(network)) / 2
  }

  if (isTRUE(scaleScores)) {
    max_abs_value <- max(abs(network), na.rm = TRUE)
    if (is.finite(max_abs_value) && max_abs_value > 0) {
      network <- network / max_abs_value
    }
  }

  filter_edges <- q > 0 && q < 1
  if (filter_edges) {
    network <- .pcNetFilterEdges(network, q, prior_matrix)
  }

  if (verbose) {
    cli::cli_alert_info("Building sparse network matrix...")
  }

  network <- .pcNetAsSparse(network, gene_names[used], drop_zeros = filter_edges)
  if (n_used < n_genes) {
    network <- .pcNetAddEmptyGenes(network, used, gene_names)
  }

  if (verbose) {
    cli::cli_alert_success("Done!")
    cli::cli_inform("Network dimensions: {nrow(network)} x {ncol(network)}")
    cli::cli_inform("Non-zero edges: {length(network@x)}")
  }

  network
}


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Sparse genes x genes matrix of prior edges. Edges listed more than once are
# counted more than once (their entry is 2, 3, ...), as in the original
# implementation.
.pcNetPriorMatrix <- function(priorNetwork, gene_names) {
  if (ncol(priorNetwork) != 2) {
    cli::cli_abort("Prior network needs to be a two column data.frame with regulators and targets")
  }
  known <- priorNetwork[[1]] %in% gene_names & priorNetwork[[2]] %in% gene_names
  priorNetwork <- priorNetwork[known, ]
  n_genes <- length(gene_names)
  Matrix::sparseMatrix(
    i = factor(priorNetwork[[1]], gene_names),
    j = factor(priorNetwork[[2]], gene_names),
    x = 1,
    dims = c(n_genes, n_genes),
    dimnames = list(gene_names, gene_names)
  )
}

# Standardizes the non-constant genes. Returns
#   X_std:        dense samples x genes matrix of the non-constant genes, each
#                 with mean 0 and standard deviation 1 (as scale(t(X))).
#   non_constant: logical, one per input gene; FALSE for genes that have the
#                 same value in every sample (their standard deviation is 0).
.pcNetStandardize <- function(X) {
  X_std <- as.matrix(Matrix::t(X))
  storage.mode(X_std) <- "double"
  dimnames(X_std) <- NULL
  n_samples <- nrow(X_std)

  # Exact test: a gene is constant if every value equals its first value
  non_constant <- colSums(X_std != rep(X_std[1L, ], each = n_samples)) > 0
  if (!all(non_constant)) {
    X_std <- X_std[, non_constant, drop = FALSE]
  }

  gene_means <- colMeans(X_std)
  X_std <- X_std - rep(gene_means, each = n_samples)

  gene_sds <- sqrt(colSums(X_std * X_std) / max(1L, n_samples - 1L))
  list(X_std = X_std / rep(gene_sds, each = n_samples),
       non_constant = non_constant)
}

# Eigendecomposition of the Gram matrix X_std X_std' = Q diag(d) Q'.
# Returns
#   d: eigenvalues in decreasing order (squared singular values of X_std)
#   Z: r x genes matrix, Z = Q' X_std. Column k holds gene k in the eigenbasis.
# The smaller of the two Gram matrices is decomposed; both give the same d
# and Z.
.pcNetEigenBasis <- function(X_std) {
  if (nrow(X_std) <= ncol(X_std)) {
    # Fewer samples than genes: samples x samples Gram matrix
    eig <- eigen(tcrossprod(X_std), symmetric = TRUE)
    d <- pmax(eig$values, 0)
    Z <- crossprod(eig$vectors, X_std)
  } else {
    # Fewer genes than samples: genes x genes matrix, then Z = diag(sqrt(d)) V'
    eig <- eigen(crossprod(X_std), symmetric = TRUE)
    d <- pmax(eig$values, 0)
    Z <- sqrt(d) * t(eig$vectors)
  }
  list(Z = Z, d = d)
}

# Dense genes x genes matrix of regression coefficients (diagonal not yet
# cleared). Row k holds the coefficients of the regression of gene k.
# Genes are independent, so with nCores > 1 they are split into blocks and
# each block of rows is computed in its own process.
.pcNetCoefficients <- function(Z, d, nComp, nCores) {
  coefficient_rows <- function(genes, Z, d, nComp, secular_weights) {
    crossprod(secular_weights(Z[, genes, drop = FALSE], d, nComp), Z)
  }

  n_genes <- ncol(Z)
  if (nCores == 1L) {
    return(coefficient_rows(seq_len(n_genes), Z, d, nComp, .pcNetSecularWeights))
  }

  # Workers receive every function they need as an argument; detaching the
  # functions from this environment keeps the data sent to workers small.
  secular_weights <- .pcNetSecularWeights
  environment(secular_weights) <- baseenv()
  environment(coefficient_rows) <- baseenv()
  gene_blocks <- split(seq_len(n_genes),
                       cut(seq_len(n_genes), nCores, labels = FALSE))

  if (.Platform$OS.type == "windows") {
    cluster <- parallel::makePSOCKcluster(nCores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    row_blocks <- parallel::parLapply(
      cluster, gene_blocks, coefficient_rows,
      Z = Z, d = d, nComp = nComp, secular_weights = secular_weights
    )
  } else {
    row_blocks <- parallel::mclapply(
      gene_blocks, coefficient_rows,
      Z = Z, d = d, nComp = nComp, secular_weights = secular_weights,
      mc.cores = nCores, mc.preschedule = TRUE
    )
    # mclapply returns an error object instead of a matrix for failed workers
    failed <- !vapply(row_blocks, is.matrix, logical(1))
    if (any(failed)) {
      messages <- unique(unlist(lapply(row_blocks[failed], as.character)))
      stop("Parallel computation failed: ", paste(messages, collapse = "; "))
    }
  }

  do.call(rbind, row_blocks)
}

# Sets edges whose absolute weight is below the q-th quantile to zero. Edges
# in the prior network keep their weight, multiplied by the number of times
# they are listed in the prior.
.pcNetFilterEdges <- function(network, q, prior_matrix = NULL) {
  abs_network <- abs(network)
  threshold <- stats::quantile(abs_network, q, na.rm = TRUE)

  # Prior edge weights, taken before filtering. Positions are linear indices
  # into the dense matrix, read from the column-compressed prior matrix.
  if (!is.null(prior_matrix)) {
    n_genes <- nrow(network)
    prior_rows <- prior_matrix@i + 1
    prior_cols <- rep.int(seq_len(n_genes) - 1, diff(prior_matrix@p))
    prior_index <- prior_rows + prior_cols * n_genes
    prior_weight <- prior_matrix@x * network[prior_index]
  }

  network[abs_network < threshold] <- 0

  if (!is.null(prior_matrix)) {
    restore <- prior_weight != 0
    network[prior_index[restore]] <- prior_weight[restore]
  }

  network
}

# Converts the dense network (zero diagonal) to a dgCMatrix. The layout
# matches the original implementation:
#   drop_zeros = FALSE: every off-diagonal entry is stored, even if zero, and
#                       the diagonal is not stored.
#   drop_zeros = TRUE:  only nonzero entries are stored (used after filtering).
.pcNetAsSparse <- function(network, gene_names, drop_zeros) {
  n_genes <- nrow(network)
  dim_names <- if (is.null(gene_names)) list(NULL, NULL) else list(gene_names, gene_names)

  if (drop_zeros) {
    # generalMatrix keeps a dgCMatrix even when the network is symmetric
    network <- methods::as(methods::as(network, "CsparseMatrix"), "generalMatrix")
    network@Dimnames <- dim_names
    return(network)
  }

  # Linear positions of the diagonal, removed from the column-major values
  diagonal <- seq.int(1L, by = n_genes + 1L, length.out = n_genes)
  methods::new(
    "dgCMatrix",
    i = rep.int(seq_len(n_genes) - 1L, n_genes)[-diagonal],
    p = seq.int(0L, by = n_genes - 1L, length.out = n_genes + 1L),
    x = as.vector(network)[-diagonal],
    Dim = c(n_genes, n_genes),
    Dimnames = dim_names
  )
}

# Expands the network of the non-constant genes to all genes. Constant genes
# get empty rows and columns (no stored entries).
.pcNetAddEmptyGenes <- function(network, used, gene_names) {
  n_genes <- length(used)
  used_index <- which(used)
  entries_per_column <- integer(n_genes)
  entries_per_column[used_index] <- diff(network@p)
  dim_names <- if (is.null(gene_names)) list(NULL, NULL) else list(gene_names, gene_names)
  methods::new(
    "dgCMatrix",
    i = used_index[network@i + 1L] - 1L,
    p = c(0L, cumsum(entries_per_column)),
    x = network@x,
    Dim = c(n_genes, n_genes),
    Dimnames = dim_names
  )
}

# Solves the per-gene eigenproblems through the secular equation.
#
# Input
#   Z:     r x genes. Column k is z = Q' y for gene k (see .pcNetEigenBasis).
#   d:     the r eigenvalues of the Gram matrix, in decreasing order.
#   nComp: number of principal components.
#
# Output
#   W: r x genes. Column k is w = sum_i u_i (u_i' z) / mu_i, where (mu_i, u_i)
#      are the top nComp eigenpairs of D - z z'. The coefficients of gene k
#      are then t(Z) %*% w.
#
# Method
#   Removing gene k lowers every eigenvalue, and the i-th largest new
#   eigenvalue mu_i lies between d[i + 1] and d[i]. It is the only root in that
#   interval of
#
#     f(mu) = 1 - sum_j z_j^2 / (d_j - mu),
#
#   which is decreasing there. Its eigenvector is proportional to
#   v = z / (d - mu), and since f(mu) = 0 gives u_i' z = 1 / ||v||, each term
#   of w is v / (mu * ||v||^2).
#
#   Each root is found as in LAPACK's dlaed4:
#   * mu is written as origin + tau, where origin is the closer of the two
#     interval ends. Then d_j - mu = (d_j - origin) - tau is computed without
#     cancellation, even when mu is extremely close to an eigenvalue.
#   * Each step replaces the sums over eigenvalues above and below the
#     interval by simple one-pole functions that match their value and slope,
#     and solves the resulting quadratic. This converges in a few steps.
#   * A bracket around the root is kept, and a bisection step is used whenever
#     the quadratic step would leave it, so the iteration always converges.
#   All genes are solved together, as vectorized column operations; genes are
#   dropped from the working set once converged.
.pcNetSecularWeights <- function(Z, d, nComp) {
  eps <- .Machine$double.eps
  n_eig_input <- nrow(Z)
  n_eig <- n_eig_input
  n_genes <- ncol(Z)

  # With very few samples, every interval needs a lower end: the remaining
  # eigenvalues of the Gram matrix are zero.
  if (n_eig < nComp + 1L) {
    n_pad <- nComp + 1L - n_eig
    Z <- rbind(Z, matrix(0, n_pad, n_genes))
    d <- c(d, numeric(n_pad))
    n_eig <- nComp + 1L
  }

  # A zero z_j means d_j is still an eigenvalue after removing the gene, with
  # an eigenvector that contributes nothing. Rather than handle that case
  # separately, raise such entries to a negligible size, which gives the same
  # result through the general formula.
  tiny <- sqrt(max(d[1L], .Machine$double.xmin)) * 1e-60
  near_zero <- abs(Z) < tiny
  if (any(near_zero)) {
    Z[near_zero] <- ifelse(Z[near_zero] < 0, -tiny, tiny)
  }
  Z2 <- Z * Z

  W <- matrix(0, n_eig, n_genes)

  for (i in seq_len(nComp)) {
    d_upper <- d[i]
    d_lower <- d[i + 1L]
    gap <- d_upper - d_lower

    # Equal eigenvalues: mu_i equals them and its eigenvector is orthogonal to
    # z, so it contributes nothing to w.
    if (!(gap > 0)) next

    rows_above <- seq_len(i)            # eigenvalues >= d_upper
    rows_below <- (i + 1L):n_eig        # eigenvalues <= d_lower

    # Pick the closer end as origin: f(midpoint) >= 0 means the root is in the
    # upper half of the interval.
    midpoint <- d_lower + gap / 2
    near_upper <- (1 - colSums(Z2 / (d - midpoint))) >= 0
    origin <- ifelse(near_upper, d_upper, d_lower)

    # Everything below is relative to origin (tau = mu - origin), per gene
    pole_lower <- ifelse(near_upper, -gap, 0)
    pole_upper <- ifelse(near_upper, 0, gap)
    bracket_lo <- ifelse(near_upper, -gap / 2, 0)
    bracket_hi <- ifelse(near_upper, 0, gap / 2)
    tau <- (bracket_lo + bracket_hi) / 2
    d_minus_origin <- outer(d, origin, "-")

    active <- seq_len(n_genes)
    for (iteration in seq_len(300L)) {
      tau_a <- tau[active]
      gaps <- d_minus_origin[, active, drop = FALSE] - rep(tau_a, each = n_eig)
      terms <- Z2[, active, drop = FALSE] / gaps     # z_j^2 / (d_j - mu)
      slopes <- terms / gaps                         # z_j^2 / (d_j - mu)^2

      sum_above <- colSums(terms[rows_above, , drop = FALSE])
      slope_above <- colSums(slopes[rows_above, , drop = FALSE])
      sum_below <- colSums(terms[rows_below, , drop = FALSE])
      slope_below <- colSums(slopes[rows_below, , drop = FALSE])
      f <- 1 - sum_above - sum_below

      # f is decreasing: f > 0 means the root is above tau, f < 0 below
      lo <- bracket_lo[active]
      hi <- bracket_hi[active]
      lo[f > 0] <- tau_a[f > 0]
      hi[f < 0] <- tau_a[f < 0]
      bracket_lo[active] <- lo
      bracket_hi[active] <- hi

      # Model: sum_below ~ a + b / (pole_lower - x), sum_above ~ c + e / (pole_upper - x),
      # matching value and slope at tau. Solving model = 1 for x = tau + step
      # gives   A step^2 - B step + C = 0.
      to_lower <- pole_lower[active] - tau_a
      to_upper <- pole_upper[active] - tau_a
      A <- f + slope_below * to_lower + slope_above * to_upper
      B <- A * (to_lower + to_upper) -
        slope_below * to_lower * to_lower - slope_above * to_upper * to_upper
      C <- f * to_lower * to_upper

      # Numerically stable quadratic roots; keep the one between the poles
      sqrt_disc <- sqrt(pmax(B * B - 4 * A * C, 0))
      half <- (B + ifelse(B >= 0, sqrt_disc, -sqrt_disc)) / 2
      root_1 <- half / A
      root_2 <- C / half
      step <- ifelse(!is.na(root_1) & root_1 > to_lower & root_1 < to_upper,
                     root_1, root_2)

      # Fall back to bisection if the step leaves the bracket
      tau_new <- tau_a + step
      outside <- is.na(tau_new) | tau_new <= lo | tau_new >= hi
      tau_new[outside] <- (lo[outside] + hi[outside]) / 2

      exact <- f == 0
      tau_new[exact] <- tau_a[exact]
      converged <- exact |
        abs(tau_new - tau_a) <= 4 * eps * abs(tau_new) |
        (hi - lo) <= 4 * eps * pmax(abs(lo), abs(hi))

      tau[active] <- tau_new
      active <- active[!converged]
      if (length(active) == 0L) break
    }

    # Add this eigenpair's term v / (mu * ||v||^2) with v = z / (d - mu)
    mu <- origin + tau
    v <- Z / (d_minus_origin - rep(tau, each = n_eig))
    W <- W + v * rep(1 / (mu * colSums(v * v)), each = n_eig)
  }

  W[seq_len(n_eig_input), , drop = FALSE]
}
