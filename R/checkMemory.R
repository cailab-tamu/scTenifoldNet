#' @export checkMemory
#' @importFrom ps ps_system_memory
#' @importFrom cli cli_warn
#' @title Checks whether there is enough memory for a number of genes
#' @description Estimates the peak memory that network construction, tensor
#'   decomposition and manifold alignment will need for a given number of
#'   genes, and compares it with the memory currently available. When the
#'   estimate does not fit, a warning reports the largest number of genes that
#'   fits in the available memory. \code{scTenifoldNet()} and
#'   \code{scTenifoldKnk()} call it after quality control, before building the
#'   networks.
#' @param nGenes An integer value. The number of genes in the networks.
#' @param nNet An integer value. The number of networks. Default: 10.
#' @param nConditions An integer value. The number of conditions compared: 2 for
#'   \code{scTenifoldNet()}, which builds and decomposes one tensor per
#'   condition, and 1 for \code{scTenifoldKnk()}. Default: 1.
#' @param warn A boolean value (TRUE/FALSE). If TRUE (default), a warning is
#'   raised when the estimated peak memory exceeds the available memory.
#' @return Invisibly, a list with \code{required} (estimated peak memory in
#'   bytes), \code{available} (memory currently available in bytes, as reported
#'   by \code{ps::ps_system_memory()}; \code{NA} if it cannot be determined),
#'   \code{total} (total memory in bytes) and \code{maxGenes} (the largest
#'   number of genes whose estimate fits in the available memory).
#' @details Memory grows with the square of the number of genes: the tensor of
#'   networks has \code{nGenes x nGenes x nNet} entries, and the decomposition
#'   keeps working copies of it. The estimate is calibrated from benchmarks of
#'   the full pipelines and is approximate; available memory also changes while
#'   R runs, so treat the result as a guide.
#' @examples
#' # Estimated peak memory for 5,000 genes and 10 networks
#' mem <- checkMemory(nGenes = 5000, nNet = 10, warn = FALSE)
#' mem$required / 1e9
#'
#' # Largest number of genes that fits in the memory available now
#' mem$maxGenes
checkMemory <- function(nGenes, nNet = 10, nConditions = 1, warn = TRUE) {
  if (!nConditions %in% c(1, 2)) {
    stop("'nConditions' must be 1 or 2")
  }
  model <- .memoryModel[[nConditions]]
  perGene2 <- model$perNetwork * nNet + model$perGene2
  required <- model$base + perGene2 * nGenes^2

  systemMemory <- tryCatch(ps::ps_system_memory(), error = function(e) NULL)
  available <- if (is.null(systemMemory)) NA_real_ else as.numeric(systemMemory$avail)
  total <- if (is.null(systemMemory)) NA_real_ else as.numeric(systemMemory$total)
  maxGenes <- if (is.na(available)) NA_real_ else
    floor(sqrt(max(available - model$base, 0) / perGene2))

  if (isTRUE(warn) && !is.na(available) && required > available) {
    gb <- function(x) format(round(x / 1e9, 1), nsmall = 1)
    cli::cli_warn(c(
      "!" = paste0("Estimated peak memory for {format(nGenes, big.mark = ',')} genes and ",
                   "{nNet} network{?s} is ~{gb(required)} GB, but only {gb(available)} GB ",
                   "of {gb(total)} GB is available."),
      "i" = paste0("With the memory available, at most ~{format(maxGenes, big.mark = ',')} ",
                   "genes can be used. Consider reducing the number of genes ",
                   "(e.g. with {.arg qc_minPCT}) or of networks ({.arg nc_nNet}).")
    ))
  }

  invisible(list(required = required, available = available, total = total,
                 maxGenes = maxGenes))
}

# Peak memory model, in bytes: base + (perNetwork * nNet + perGene2) * nGenes^2.
# Fitted to the peak R heap (gc() "max used") of full pipeline runs with
# 1,000-5,000 genes, 5-20 networks and 300-5,000 cells; predictions are within
# 5% (one condition) and 10% (two conditions) above 2,000 genes. The number of
# cells has little effect because each network uses a fixed subsample.
.memoryModel <- list(
  # scTenifoldKnk: one tensor
  list(base = 2.0e8, perNetwork = 29.3, perGene2 = 43.5),
  # scTenifoldNet: two tensors
  list(base = 3.3e8, perNetwork = 37.6, perGene2 = 82.2)
)
