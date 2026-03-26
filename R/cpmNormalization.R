#' @export cpmNormalization
#' @title Performs counts per million (CPM) data normalization
#' @importFrom cli cli_alert_info cli_alert_success
#' @description This function normalizes the count data present in a given
#'   matrix using counts per million normalization (CPM). Each gene count for
#'   each cell is divided by the total counts for that cell and multiplied by
#'   1e6. No log-transformation is applied.
#' @param X Raw counts matrix with cells as columns and genes (symbols) as rows.
#' @param label Optional character label (e.g. "X", "Y") prepended to
#'   progress messages when running inside a pipeline.
#' @return A dgCMatrix object with the count per million (CPM) normalized
#'   values.
#' @references Vallejos, Catalina A., et al. "Normalizing single-cell RNA
#'   sequencing data: challenges and opportunities." Nature methods 14.6 (2017):
#'   565.
#' @examples 
#' library(scTenifoldNet)
#' 
#' # Simulating of a dataset following a negative binomial distribution with high sparcity (~67%)
#' nCells = 2000
#' nGenes = 100
#' set.seed(1)
#' X <- rnbinom(n = nGenes * nCells, size = 20, prob = 0.98)
#' X <- round(X)
#' X <- matrix(X, ncol = nCells)
#' rownames(X) <- c(paste0('ng', 1:90), paste0('mt-', 1:10))
#' 
#' # Performing Single cell quality control
#' qcOutput <- scQC(
#'   X = X,
#'   minLibSize = 30,
#'   removeOutlierCells = TRUE,
#'   minPCT = 0.05,
#'   maxMTratio = 0.1
#' )
#' # Performing Counts per million Normalization (CPM)
#' normalizationOutput <- cpmNormalization(qcOutput)
#' 
#' # Visualizing the differences
#' oldPar <- par(no.readonly = TRUE)
#' 
#' par(
#'   mfrow = c(1, 2),
#'   mar = c(3, 3, 1, 1),
#'   mgp = c(1.5, 0.5, 0)
#' )
#' plot(
#'   Matrix::colSums(qcOutput),
#'   ylab = 'Library Size',
#'   xlab = 'Cell',
#'   main = 'Before CPM Normalization'
#' )
#' plot(
#'   Matrix::colSums(normalizationOutput),
#'   ylab = 'Library Size',
#'   xlab = 'Cell',
#'   main = 'After CPM Normalization'
#' )
#' 
#' par(oldPar)

cpmNormalization <- function(X, label = NULL) {
  nGenes <- nrow(X)
  nCells <- ncol(X)
  tag <- if (!is.null(label)) paste0("[", label, "] ") else ""
  cli::cli_alert_info("{tag}CPM normalization: {nGenes} genes x {nCells} cells")
  X <- Matrix::t(Matrix::t(X) / Matrix::colSums(X)) * 1e6
  X <- as(X, 'dgCMatrix')
  cli::cli_alert_success("{tag}CPM normalization complete")
  return(X)
}

