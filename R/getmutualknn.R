#'
#' @title getmutualknn
#' @description Get a mutual kNN graph (M-kNN) from the latent space (LSI) of single cell profiles
#'
#' @param lsimat A dense or sparse matrix of cell (rows) x LSI (cols)
#' @param num_k A numeric input (30 by default) for k used for graph construction
#'
#' @return A sparse matrix indicating adjacent matrix from single cells used for graph construction
#' @export
#'
#' @import igraph
#' @import RANN
#' @importFrom dplyr %>%
#'
#' @examples
#' \dontrun{
#' mknn_graph <- getmutualknn(lsimat=LSI_mat, num_k=30)
#' }
#'
#'
getmutualknn <- function(lsimat=lsimat, num_k=30){

  stopifnot("num_k must be numeric" = is.numeric(num_k))
  mymat=lsimat # this is important
  K=num_k
  # calculate euclidean distances between cells
  message("[info] fast knn")
  knn.info<- RANN::nn2(mymat, k=K)
  ## convert to adjacancy matrix
  knn <- knn.info$nn.idx
  tempx <- rep(0, nrow(knn))
  for (i in 1:nrow(knn)){
    xx <- apply(knn[knn[i, -1], ],
                1,
                function(x) {any(x==i)}
    )
    tempx[i] <- sum(xx)
  }

  knn2 <- list()
  length(knn2) <- nrow(knn)
  for (i in 1:nrow(knn)){
    xx <- apply(knn[knn[i, -1], ],
                1,
                function(x) {any(x==i)}
    )
    if(sum(xx)>0){
      temp_knn <- knn[i, c(TRUE, xx)]
      temp_el <- cbind(temp_knn[1], c(temp_knn[-1]))
    } else {
      temp_el <- knn[i, 1:2]
    }
    knn2[[i]] <- temp_el
  }
  el <- do.call(rbind.data.frame, knn2) %>% as.matrix

  adj <- igraph::get.adjacency(igraph::graph.edgelist(el))
  mutualknn <- 1*((adj + t(adj)) > 0)
  colnames(mutualknn) <- rownames(mutualknn) <- rownames(lsimat)
  return(mutualknn)
}
