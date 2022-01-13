#' Calculate AUC as Mann–Whitney U statistic from a vector of ranks
#' 
#' @param rank_value    the ranks
#' @param maxRank       limit of ranks to include 
#' @param sparse        is the vector sparse
#' 
#' @return              AUC value for the vector
u_stat <- function(rank_value, maxRank=1000, sparse=F){

  if(sparse==T){
    rank_value[rank_value==0] <- maxRank+1
  }

  insig <- rank_value > maxRank
  if(all(insig)) {
    return(0L)
  } else {
    rank_value[insig] <- maxRank+1
    rank_sum = sum(rank_value)
    len_sig <- length(rank_value)
    u_value = rank_sum - (len_sig * (len_sig + 1))/2
    auc = 1 - u_value/(len_sig * maxRank)
    return(auc)
  }
}

