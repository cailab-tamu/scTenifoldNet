#' @title Canonical Polyadic Decomposition
#'
#' @description Canonical Polyadic (CP) decomposition of a tensor, aka
#'   CANDECOMP/PARAFAC. Approximate a K-Tensor using a sum of
#'   \code{num_components} rank-1 K-Tensors. A rank-1 K-Tensor can be written
#'   as an outer product of K vectors. There are a total of
#'   \code{num_components * tnsr$num_modes} vectors in the output, stored in
#'   \code{tnsr$num_modes} matrices, each with \code{num_components} columns.
#'   This is an iterative algorithm, with two possible stopping conditions:
#'   either relative error in Frobenius norm has gotten below \code{tol}, or
#'   the \code{max_iter} number of iterations has been reached. For more details
#'   on CP decomposition, consult Kolda and Bader (2009).
#' @details Uses the Alternating Least Squares (ALS) estimation procedure. A
#'   progress bar is included to help monitor operations on large tensors.
#' @name cpDecomposition
#' @importFrom methods is
#' @importFrom stats rnorm
#' @importFrom utils tail
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @param tnsr Tensor with K modes.
#' @param num_components The number of rank-1 K-Tensors to use in approximation.
#' @param max_iter Maximum number of iterations if error stays above \code{tol}.
#' @param tol Relative Frobenius norm error tolerance.
#' @return A list containing the following \describe{
#'   \item{\code{lambdas}}{A vector of normalizing constants, one for each component.}
#'   \item{\code{U}}{A list of matrices, one for each mode, each with \code{num_components} columns.}
#'   \item{\code{conv}}{Whether or not \code{resid} < \code{tol} by the last iteration.}
#'   \item{\code{norm_percent}}{The percent of Frobenius norm explained by the approximation.}
#'   \item{\code{est}}{Estimate of \code{tnsr} after compression.}
#'   \item{\code{fnorm_resid}}{The Frobenius norm of the error.}
#'   \item{\code{all_resids}}{Vector containing the Frobenius norm of error for all iterations.}
#' }
#' @references T. Kolda, B. Bader, "Tensor decomposition and applications".
#'   SIAM Applied Mathematics and Applications 2009.

cpDecomposition <- function(tnsr, num_components = NULL, max_iter = 25,
                            tol = 1e-5) {
  if (is.null(num_components)) stop("num_components must be specified")
  stopifnot(is(tnsr, "list"))

  num_modes <- tnsr$num_modes
  modes <- tnsr$modes
  tnsr_data <- tnsr$data
  tnsr_norm_sq <- sum(tnsr_data * tnsr_data)
  tnsr_norm <- sqrt(tnsr_norm_sq)
  R <- num_components

  # Initialize factor matrices randomly
  U_list <- vector("list", num_modes)
  for (m in seq_len(num_modes)) {
    U_list[[m]] <- matrix(rnorm(modes[m] * R), nrow = modes[m], ncol = R)
  }

  curr_iter <- 1
  converged <- FALSE
  fnorm_resid <- rep(0, max_iter)
  lambdas <- numeric(R)
  prev_resid <- Inf

  pb <- cli::cli_progress_bar("CP decomposition", total = max_iter)

  if (num_modes == 3L) {
    # ================================================================
    # FAST PATH: 3-mode tensor (the only case used in scTenifoldNet)
    # Uses slicewise MTTKRP (avoids forming large Khatri-Rao products)
    # and algebraic convergence check (avoids full tensor reconstruction)
    # ================================================================
    I <- modes[1]; J <- modes[2]; K <- modes[3]

    # Pre-extract slices for cache-friendly access
    slices <- vector("list", K)
    for (k in seq_len(K)) slices[[k]] <- tnsr_data[,, k]

    while ((curr_iter < max_iter) && (!converged)) {
      cli::cli_progress_update(id = pb, set = curr_iter)

      # --- Mode 1 update ---
      V <- crossprod(U_list[[2]]) * crossprod(U_list[[3]])
      mttkrp <- matrix(0, I, R)
      for (k in seq_len(K))
        mttkrp <- mttkrp + slices[[k]] %*% U_list[[2]] * rep(U_list[[3]][k,], each = I)
      tmp <- mttkrp %*% solve(V)
      lambdas <- colSums(abs(tmp))
      U_list[[1]] <- sweep(tmp, 2, lambdas, "/")

      # --- Mode 2 update ---
      V <- crossprod(U_list[[1]]) * crossprod(U_list[[3]])
      mttkrp <- matrix(0, J, R)
      for (k in seq_len(K))
        mttkrp <- mttkrp + crossprod(slices[[k]], U_list[[1]]) * rep(U_list[[3]][k,], each = J)
      tmp <- mttkrp %*% solve(V)
      lambdas <- colSums(abs(tmp))
      U_list[[2]] <- sweep(tmp, 2, lambdas, "/")

      # --- Mode 3 update ---
      V <- crossprod(U_list[[1]]) * crossprod(U_list[[2]])
      mttkrp <- matrix(0, K, R)
      for (k in seq_len(K))
        mttkrp[k,] <- colSums(U_list[[1]] * (slices[[k]] %*% U_list[[2]]))
      tmp <- mttkrp %*% solve(V)
      lambdas <- colSums(abs(tmp))
      U_list[[3]] <- sweep(tmp, 2, lambdas, "/")

      # Algebraic convergence: ||X-est||^2 = ||X||^2 - 2<X,est> + ||est||^2
      inner <- sum(lambdas * colSums(U_list[[3]] * mttkrp))
      Gamma <- crossprod(U_list[[1]]) * crossprod(U_list[[2]]) * crossprod(U_list[[3]])
      est_norm_sq <- sum(outer(lambdas, lambdas) * Gamma)
      curr_resid <- sqrt(max(tnsr_norm_sq - 2 * inner + est_norm_sq, 0))

      fnorm_resid[curr_iter] <- curr_resid
      if (curr_iter > 1 && abs(curr_resid - prev_resid) / tnsr_norm < tol) {
        converged <- TRUE
        cli::cli_progress_update(id = pb, set = max_iter, force = TRUE)
      } else {
        prev_resid <- curr_resid
        curr_iter <- curr_iter + 1
      }
    }

    if (!converged) cli::cli_progress_update(id = pb, set = max_iter, force = TRUE)
    cli::cli_progress_done(id = pb)

    # Reconstruct tensor ONCE at the end
    est_data <- array(0, dim = modes)
    for (k in seq_len(K))
      est_data[,, k] <- tcrossprod(
        sweep(U_list[[1]], 2, lambdas * U_list[[3]][k, ], "*"),
        U_list[[2]]
      )

  } else {
    # ================================================================
    # GENERAL PATH: N-mode tensor (fallback)
    # ================================================================

    # Pre-compute mode-m unfoldings
    unfolded_mat <- vector("list", num_modes)
    for (m in seq_len(num_modes)) {
      perm <- c(m, seq_len(num_modes)[-m])
      mat <- aperm(tnsr_data, perm)
      dim(mat) <- c(modes[m], prod(modes[-m]))
      unfolded_mat[[m]] <- mat
    }

    khatri_rao_rev <- function(L) {
      L <- rev(L)
      nc <- ncol(L[[1]])
      nr <- vapply(L, nrow, integer(1))
      retmat <- matrix(0, nrow = prod(nr), ncol = nc)
      for (j in seq_len(nc)) {
        col_j <- L[[1]][, j]
        for (i in 2:length(L)) col_j <- kronecker(col_j, L[[i]][, j])
        retmat[, j] <- col_j
      }
      retmat
    }

    while ((curr_iter < max_iter) && (!converged)) {
      cli::cli_progress_update(id = pb, set = curr_iter)
      for (m in seq_len(num_modes)) {
        V <- Reduce("*", lapply(U_list[-m], crossprod))
        KR <- khatri_rao_rev(U_list[-m])
        mttkrp <- unfolded_mat[[m]] %*% KR
        tmp <- mttkrp %*% solve(V)
        lambdas <- colSums(abs(tmp))
        U_list[[m]] <- sweep(tmp, 2, lambdas, "/")
      }

      # Algebraic convergence
      inner <- sum(lambdas * colSums(U_list[[num_modes]] * mttkrp))
      Gamma <- Reduce("*", lapply(U_list, crossprod))
      est_norm_sq <- sum(outer(lambdas, lambdas) * Gamma)
      curr_resid <- sqrt(max(tnsr_norm_sq - 2 * inner + est_norm_sq, 0))

      fnorm_resid[curr_iter] <- curr_resid
      if (curr_iter > 1 && abs(curr_resid - prev_resid) / tnsr_norm < tol) {
        converged <- TRUE
        cli::cli_progress_update(id = pb, set = max_iter, force = TRUE)
      } else {
        prev_resid <- curr_resid
        curr_iter <- curr_iter + 1
      }
    }

    if (!converged) cli::cli_progress_update(id = pb, set = max_iter, force = TRUE)
    cli::cli_progress_done(id = pb)

    KR <- khatri_rao_rev(U_list[-1])
    est_data <- array(U_list[[1]] %*% (lambdas * t(KR)), dim = modes)
  }

  fnorm_resid <- fnorm_resid[fnorm_resid != 0]
  norm_percent <- (1 - (tail(fnorm_resid, 1) / tnsr_norm)) * 100
  est <- list(data = est_data, modes = modes, num_modes = num_modes)

  invisible(list(lambdas = lambdas, U = U_list, conv = converged, est = est,
                 norm_percent = norm_percent, fnorm_resid = tail(fnorm_resid, 1),
                 all_resids = fnorm_resid))
}