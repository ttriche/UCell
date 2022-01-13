#' Get signature scores from precomputed rank matrix
#' 
#' @param     ranks_matrix    a rank matrix
#' @param     features        signature features
#' @param     chunk.size      how many cells per matrix chunk
#' @param     w_neg           can W be negative? 
#' @param     ncores          how many cores to use for parallelization
#' @param     force.gc        force garbage collection to recover RAM? (FALSE)
#' @param     name            name suffix of temporary objects ("_UCell")
#' 
#' @return                    a list of results 
#' 
#' @import    data.table
#' @import    future    
rankings2Uscore <- function(ranks_matrix, features, chunk.size=1000, w_neg=1,
                            ncores=1, force.gc=FALSE, name="_UCell") {

  #Check if all genes in signatures are present in the stored signatures
  ranks_matrix <- check_genes(ranks_matrix, features)

  maxRank <- max(ranks_matrix)
  split.data <- split_data.matrix(matrix=ranks_matrix, chunk.size=chunk.size)
  rm(ranks_matrix)

  if (ncores>1) {
    future::plan(future::multisession(workers=future_param_ncores))

    meta.list <- future_lapply(
      X = split.data,
      FUN = function(x) {

        dense <- as.matrix(x)
        dense <- as.data.table(dense, keep.rownames=TRUE)
        setkey(dense, rn, physical=F)

        cells_AUC <- u_stat_signature_list(features, dense, maxRank=maxRank, sparse=T, w_neg=w_neg)
        colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)
        
        if (force.gc) {
          dense <- NULL
          gc()
        }
        return(list(cells_AUC=cells_AUC))
      },
      future.seed = future_param_seed
    )
    future::plan(strategy = "sequential")

  } else {

    meta.list <- lapply(
      X = split.data,
      FUN = function(x) {

        dense <- as.matrix(x)
        dense <- as.data.table(dense, keep.rownames=TRUE)
        setkey(dense, rn, physical=F)

        cells_AUC <- u_stat_signature_list(features, dense, maxRank=maxRank, sparse=T, w_neg=w_neg)
        colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)
        
        if (force.gc) {
          dense <- NULL
          gc()
        }

        return(list(cells_AUC=cells_AUC))
      } )
  }
  return(meta.list)
}
