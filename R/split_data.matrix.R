#' split data matrix into cell chunks
#'
#' @param   matrix      the data matrix 
#' @param   chunk.size  how many cells per submatrix 
#' 
#' @return  a list of chunked matrices 
split_data.matrix <- function(matrix, chunk.size=1000) {
  ncols <- dim(matrix)[2]
  nchunks <- (ncols-1) %/% chunk.size + 1

  split.data <- list()
  min <- 1
  for (i in 1:nchunks) {
    if (i == nchunks-1) {  #make last two chunks of equal size
       left <- ncols-(i-1)*chunk.size
       max <- min+round(left/2)-1
    } else {
       max <- min(i*chunk.size, ncols)
    }
    split.data[[i]] <- matrix[,min:max]
    min <- max+1    #for next chunk
  }
  return(split.data)
}
