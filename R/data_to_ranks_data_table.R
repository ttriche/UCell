#' Calculate features' ranks from expression data matrices
#' 
#' @param data              the expression data matrix 
#' @param ties.method       how to break ties
#' 
#' @return                  a data.table of ranks 
#' 
#' @import data.table
data_to_ranks_data_table <- function(data, ties.method="average") {
  dt <- as.data.table(as.matrix(data))
  rnaDT.ranks.dt <- dt[, lapply(.SD, function(x) frankv(x,ties.method=ties.method,order=c(-1L)))]
  rnaDT.ranks.rownames <- rownames(data)
  rnaDT.ranks.dt.rn <- cbind(rn=rnaDT.ranks.rownames, rnaDT.ranks.dt)
  setkey(rnaDT.ranks.dt.rn, rn, physical = F)
  return(rnaDT.ranks.dt.rn)
}


