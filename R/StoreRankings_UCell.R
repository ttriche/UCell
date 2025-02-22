#' Calculate and store gene rankings for a single-cell dataset
#'
#' Given a gene vs. cell matrix, calculates the rankings of expression for all genes in each cell. 
#' 
#' While \code{\link{ScoreSignatures_UCell}} can be used 'on the fly' to evaluate signatures in a query dataset, it requires recalculating gene
#' ranks at every execution. If you have a large dataset and plan to experiment with multiple signatures, evaluating the same dataset multiple times,
#' this function allows you to store pre-calculated ranks so they do not have to be recomputed every time. Pre-calculated ranks can then be applied to the
#' function \code{\link{ScoreSignatures_UCell}} to evaluate gene signatures in a significantly faster way on successive iterations.
#'
#' @param matrix A gene vs. cell data matrix, either in sparse or dense format
#' @param maxRank Maximum number of genes to rank per cell; above this rank, a given gene is considered as not expressed
#' @param chunk.size Number of cells to be processed simultaneously (lower size requires slightly more computation but reduces memory demands)
#' @param ncores Number of processors to parallelize computation. Requires package \code{future}
#' @param ties.method How ranking ties should be resolved (passed on to [data.table::frank])
#' @param force.gc Explicitly call garbage collector to reduce memory footprint
#' @param seed Integer seed for [future.apply::future_lapply] parallel execution
#' @return Returns a sparse matrix of pre-calculated ranks that can be used multiple times to evaluate different signatures
#' @examples
#' ## Not run:
#' library(UCell)
#' my.matrix <- UCell::sample.matrix
#' ranks <- StoreRankings_UCell(my.matrix)
#' ranks[1:5,1:5]
#' gene.sets <- list( Tcell_signature = c("CD2","CD3E","CD3D"),
#'                  Myeloid_signature = c("SPI1","FCER1G","CSF1R"))
#' scores <- ScoreSignatures_UCell(features=gene.sets, precalc.ranks=ranks)
#' ## End (Not run)
#' @export
StoreRankings_UCell <- function(matrix, maxRank=1500, chunk.size=1000, ncores=1, 
                                ties.method="average", force.gc=FALSE, seed=123) {
  
  if (ncores>1) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      message("Warning: package 'future.apply' not installed. Running UCell on single core.")
      ncores <- 1
    } else {
      require(future.apply)
      future_param_seed <<- seed
      future_param_ncores <<- ncores
    }
  }
  
  features <- rownames(matrix)[1]  #dummy signature
  meta.list <- calculate_Uscore(matrix, features=features, maxRank=maxRank, chunk.size=chunk.size,
                                ncores=ncores, ties.method=ties.method, storeRanks=T, force.gc=force.gc)
  
  ranks.all <- lapply(meta.list,function(x) rbind(x[["cells_rankings"]]))
  ranks.all <- Reduce(cbind, ranks.all)
  
  return(ranks.all)

}
