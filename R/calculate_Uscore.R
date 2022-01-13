#' Calculate rankings and scores for query data and given signature set
#' 
#' @param   matrix        the data matrix 
#' @param   features      signature features 
#' @param   maxRank       rank cutoff (1500) 
#' @param   chunk.size    cells per matrix chunk (1000) 
#' @param   ncores        number of cores to use for parallelization (1) 
#' @param   w_neg         can W be negative? (1)
#' @param   ties.method   how to break ties? ("average")
#' @param   storeRanks    store ranks? (FALSE) 
#' @param   force.gc      force garbage collection? (FALSE) 
#' @param   name          suffix for temp objects ("_UCell") 
#' 
#' @return  a list of results 
#' 
#' @import  future
#' @import  Matrix
calculate_Uscore <- function(matrix, features,  maxRank=1500, chunk.size=1000, ncores=1, w_neg=1,
                             ties.method="average", storeRanks=FALSE, force.gc=FALSE, name="_UCell") {

  #Make sure we have a sparse matrix
  if (class(matrix) != "dgCMatrix") {
    matrix <- Matrix::Matrix(as.matrix(matrix),sparse = T)
  }
  #Check if all genes in signatures are present in the data matrix
  matrix <- check_genes(matrix, features)

  #Split into manageable chunks
  split.data <- split_data.matrix(matrix=matrix, chunk.size=chunk.size)

  #Parallelize?
  if (ncores>1) {
    future::plan(future::multisession(workers=future_param_ncores))

    meta.list <- future_lapply(
      X = split.data,
      FUN = function(x) {

        cells_rankings <- data_to_ranks_data_table(x, ties.method = ties.method)
        cells_AUC <- u_stat_signature_list(features, cells_rankings, maxRank=maxRank, sparse=F, w_neg=w_neg)

        colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)

        if (storeRanks==T){
          gene.names <- as.character(as.matrix(cells_rankings[,1]))
          #make sparse
          cells_rankings[cells_rankings>maxRank] <- 0
          ranks.sparse <- Matrix::Matrix(as.matrix(cells_rankings[,-1]),sparse = T)
          dimnames(ranks.sparse)[[1]] <- gene.names
          
          if (force.gc) {
             cells_rankings <- NULL
             gc()
          }
          return(list(cells_rankings=ranks.sparse, cells_AUC=cells_AUC))

        } else {
          if (force.gc) {
            cells_rankings <- NULL
            gc()
          }
          return(list(cells_AUC=cells_AUC))
        }
      },
      future.seed = future_param_seed
    )
    future::plan(strategy = "sequential")

  } else {
    meta.list <- lapply(
      X = split.data,
      FUN = function(x) {
        cells_rankings <- data_to_ranks_data_table(x, ties.method = ties.method)
        cells_AUC <- u_stat_signature_list(features, cells_rankings, maxRank=maxRank, sparse=F, w_neg=w_neg)
        colnames(cells_AUC) <- paste0(colnames(cells_AUC),name)

        if (storeRanks==T){
          gene.names <- as.character(as.matrix(cells_rankings[,1]))
          #make sparse
          cells_rankings[cells_rankings>maxRank] <- 0
          ranks.sparse <- Matrix::Matrix(as.matrix(cells_rankings[,-1]),sparse = T)
          dimnames(ranks.sparse)[[1]] <- gene.names
          if (force.gc) {
            cells_rankings <- NULL
            gc()
          }
          return(list(cells_rankings=ranks.sparse, cells_AUC=cells_AUC))

        } else {
          if (force.gc) {
            cells_rankings <- NULL
            gc()
          }
          return(list(cells_AUC=cells_AUC))
        }

      } )
  }
  return(meta.list)
}
