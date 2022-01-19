#' @title randomWalk_sparse
#'
#' @description Calculate the network propagation score using a set of seed cells and cell-to-cell graph
#'
#' @param intM a sparse matrix indicating the adjacent matrix (m x m, where m is the cell number) of cell-to-cell network (M-kNN graph)
#' @param queryGenes a logical vector indicating seed cells (TRUE) and non-seed cells (FALSE) with length of m, where m is the cell number. The length and position are corresponding to intM
#' @param gamma a numeric value indicating the probability of node restart at each step of random walk
#' @param seedWeight a numeric vector indicating the weight assigned to each node. "NO" (by default) means considering weights are equal (no weight)
#' @param stationary_cutoff delta used for determine the stationary state between any two adjacent iterations
#'
#' @return a numeric vector of network propagation score with length of m, where m is the cell number
#' @export
#'
#' @import Matrix
#'
#' @examples
#' \dontrun{
#' randomWalk_sparse <- tfidf(intM, queryGenes, gamma=0.05, seedWeight="NO", stationary_cutoff=1e-5)}
#'
randomWalk_sparse <- function(intM, queryGenes, gamma=0.05, seedWeight="NO", stationary_cutoff=1e-5) {
       if(sum(!queryGenes %in% row.names(intM))>0) {
           stop("queryGenes contains genes not found in intMat")
       }
  Ng <- nrow(intM) # Ng is the dimension of intM

  intM <- t(t(intM) / colSums(intM))

  p0 <- numeric(length=Ng)
  names(p0) <- row.names(intM)

  # normalize the node score (sum as 1)
  p0[queryGenes] <- 1
  if(seedWeight[1]!="NO"){
    p0 <- seedWeight
    message("Will use weight for the seed")
  }
  p0 <- p0/sum(p0) # revise here if you want value the weight for each node

  rwr <- function(W, P0, gamma) {
    W <- t(W)
    PT <- P0
    k <- 0 # iteration steps
    delta <- 1 # initial the delta
    while  (delta > stationary_cutoff) {
      PT1 <- (1-gamma)*W
      PT2 <- PT1 %*% t(PT)
      PT3 <- (gamma*P0)
      PT4 <- t(PT2) + PT3
      delta <- sum(abs(PT4 - PT))
      PT <- PT4
      k <- k + 1

    }
    message("Stationary step: ", k)
    message("Stationary Delta: ", delta)
    return(drop(PT)) # drop the dgeMatrix format and return a numeric vector
  }
  res <- rwr(W=t(intM), P0=t(p0), gamma)
  return(res)
}
