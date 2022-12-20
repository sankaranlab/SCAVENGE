#' @title do_lsi
#'
#' @description Get a LSI matrix (cell-by-LSI) from TF-IDF matrix from single cell profiles
#'
#' @param mat (sparse matrix)  A saparse matrix of TF-IDF matrix (LSI) used for SVD
#' @param dims number of LSIs to calculate (default:30)
#'
#' @return a sparse matrix of LSI
#' @export
#'
#' @import irlba
#'
#' @examples
#' SE_pbmc5k <- example_data(name="pbmc5k_SE.rda")
#' peak_by_cell_mat <- SummarizedExperiment::assay(SE_pbmc5k)
#' #### Downsample to speed up example ####
#' peak_by_cell_mat <- peak_by_cell_mat[seq_len(500), seq_len(100)]
#' tfidf_mat <- tfidf(bmat=peak_by_cell_mat)
#' lsi_mat <- do_lsi(mat=tfidf_mat, dims=30)
do_lsi <- function(mat, dims=30) {
  message("SVD analysis of TF-IDF matrix")
  pca.results <- irlba::irlba(t(mat), nv=dims, fastpath=FALSE)
  PCA_result <- pca.results$u %*% diag(pca.results$d)
  rownames(PCA_result) <- colnames(mat)
  colnames(PCA_result) <- paste0('LSI_', seq_len(dims))
  return(PCA_result)
}
