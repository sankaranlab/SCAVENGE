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
#' @examples
#' \dontrun{ lsi_mat <- do_lsi(tfidf_mat, dims=30)}
#'
#' @import irlba
#'
do_lsi <- function(mat, dims=30) {
  message("SVD analysis of TF-IDF matrix")
  pca.results <- irlba::irlba(t(mat), nv=dims, fastpath=FALSE)
  PCA_result <- pca.results$u %*% diag(pca.results$d)
  rownames(PCA_result) <- colnames(mat)
  colnames(PCA_result) <- paste0('LSI_', 1:dims)
  return(PCA_result)
}
