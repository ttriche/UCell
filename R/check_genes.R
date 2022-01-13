#' Check if all genes in signatures are found in data matrix
#' otherwise add zero counts in data-matrix to complete it
#' 
#' @param matrix    the data matrix
#' @param features  features that must be present (or else completed)
#' 
#' @return          the signature-complete data matrix
check_genes <- function(matrix, features) {
   features <- unlist(features)
   features <- gsub("[-+]$","",features,perl=T)
   missing <- setdiff(features, rownames(matrix))
   ll <- length(missing)
   
   if (ll/length(features) > 0.5) {
      warning(sprintf("Over half of genes (%s%%) in specified signatures are missing from data. Check the integrity of your dataset\n", round(100*ll/length(features))))
   }
   
   if (ll>0) {
      add.mat <- Matrix::sparseMatrix(length(missing), ncol(matrix))
      rownames(add.mat) <- missing
      matrix <- rbind(matrix, add.mat)
      
      missing.concatenate <- paste(missing, collapse=",")
      warning(sprintf("The following genes were not found and will be imputed to exp=0:\n* %s",missing.concatenate))
   }
   return(matrix)
}
