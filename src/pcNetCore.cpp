// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' @title Core PC Network computation (C++ backend)
//'
//' @description Computes the gene coefficient matrix using principal component
//' regression. For each gene K, regresses it against principal components
//' derived from all other genes via truncated SVD.
//'
//' @param X_std A standardized numeric matrix (samples x genes), centered
//'   and scaled by column.
//' @param nComp Number of principal components to use.
//'
//' @return A numeric matrix (n_genes x (n_genes - 1)) of regression coefficients.
//'
//' @keywords internal
// [[Rcpp::export]]
arma::mat pcNetCoreRcpp(const arma::mat& X_std, int nComp) {
  int n_samples = X_std.n_rows;
  int n_genes   = X_std.n_cols;

  arma::mat coefficient_matrix(n_genes, n_genes - 1);

  for (int K = 0; K < n_genes; K++) {
    // Target gene to predict
    arma::vec target_gene = X_std.col(K);

    // Design matrix: all columns except K
    arma::mat design_matrix(n_samples, n_genes - 1);
    int col_idx = 0;
    for (int j = 0; j < n_genes; j++) {
      if (j != K) {
        design_matrix.col(col_idx) = X_std.col(j);
        col_idx++;
      }
    }

    // Economy SVD of design matrix, then truncate to nComp components
    arma::mat U, V;
    arma::vec s;
    arma::svd_econ(U, s, V, design_matrix);

    arma::mat principal_components = V.cols(0, nComp - 1);

    // Project design matrix onto principal components: (n_samples x nComp)
    arma::mat pc_scores = design_matrix * principal_components;

    // Normalize each PC score column by its squared norm
    arma::rowvec score_sq_norms = arma::sum(arma::square(pc_scores), 0);
    arma::mat pc_scores_normalized = pc_scores.each_row() / score_sq_norms;

    // Regress target gene on normalized PC scores (equivalent to colSums(target * normalized))
    arma::vec pc_coefficients = pc_scores_normalized.t() * target_gene;

    // Transform coefficients back to original gene space
    arma::vec gene_coefficients = principal_components * pc_coefficients;

    coefficient_matrix.row(K) = gene_coefficients.t();
  }

  return coefficient_matrix;
}
