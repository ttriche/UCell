#' Calculate AUC for a list of signatures, from a ranks matrix
#' 
#' @param   sig_list      a list of signatures
#' @param   ranks_matrix  a matrix of ranks 
#' @param   maxRank       rank cutoff, for u_stat
#' @param   sparse        is the rank matrix sparse? (FALSE) 
#' @param   w_neg         can W be na
#' 
#' @import  data.table
u_stat_signature_list <- function(sig_list, ranks_matrix, maxRank=1000, sparse=F, w_neg=1) {

  u_matrix <- sapply(sig_list, function(sig) {
     sig_neg <- grep('-$', unlist(sig), perl=T, value=T)
     sig_pos <- setdiff(unlist(sig), sig_neg)
     
     if (length(sig_pos)>0) {
        sig_pos <- gsub('\\+$','',sig_pos,perl=T)
        u_p <- as.numeric(ranks_matrix[sig_pos, lapply(.SD, function(x) u_stat(x,maxRank = maxRank,sparse=sparse)),.SDcols=-1, on="rn"])
     } else {
        u_p <- rep(0, dim(ranks_matrix)[2]-1)
     }
     if (length(sig_neg)>0) {
       sig_neg <- gsub('-$','',sig_neg,perl=T)
       u_n <- as.numeric(ranks_matrix[sig_neg, lapply(.SD, function(x) u_stat(x,maxRank = maxRank,sparse=sparse)),.SDcols=-1, on="rn"])
     } else {
       u_n <- rep(0, dim(ranks_matrix)[2]-1)
     }
     
     diff <- u_p - w_neg*u_n
     diff[diff<0] <- 0
     return(diff)
  })

  rownames(u_matrix) <- colnames(ranks_matrix)[-1]
  return (u_matrix)
}

