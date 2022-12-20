
#' @title tfidf
#'
#' @description Get a normalized TFIDF matrix from peak-by-cell matrix
#' from single cell profiles.
#'
#' @param bmat a sparse/dense matrix or dataframe indicating count matrix
#' (peak-by-cell matrix).
#' @param mat_binary a logic value to indicate if input bmat is
#' spare matrix or not. bmat will be converted into a sparse matrix.
#' @param TF a logic value to indicate if term frequency (TF) normalization
#' is performed.
#' @param log_TF a logic value to indicate if natural logarithm values of the
#'  TF are calculated and used.
#' @param scale_factor a scale factor used to multiple the resulting
#' TF-IDF matrix.
#'
#' @returns a sparse matrix of normalized TFIDF matrix.
#' @export
#'
#' @import Matrix
#' @importFrom stats quantile
#'
#' @examples
#' SE_pbmc5k <- example_data(name="pbmc5k_SE.rda")
#' peak_by_cell_mat <- SummarizedExperiment::assay(SE_pbmc5k)
#' #### Downsample to speed up example ####
#' peak_by_cell_mat <- peak_by_cell_mat[seq_len(500), seq_len(100)]
#' tfidf_mat <- tfidf(bmat=peak_by_cell_mat)
tfidf <- function(bmat,
                  mat_binary=TRUE,
                  TF=TRUE,
                  log_TF=TRUE,
                  scale_factor=100000) {
  if (mat_binary) {
    # "term frequency" method
    bmat@x[bmat@x>=1] <- 1
    message("[info] binarize matrix")
  }
  if (TF) {
    # "term frequency" method
    tf <- t(t(bmat) / Matrix::colSums(bmat))
  } else {
    # "raw count" method
    tf <- bmat
  }
  message("[info] calculate tf")
  # Either TF method can optionally be log scaled
  if (log_TF) {
    if (TF) {
      tf@x <- log1p(tf@x * scale_factor)
    } else {
      tf@x <- log1p(tf@x * 1)
    }
  }
  message("[info] calculate idf")
  # IDF "inverse document frequency smooth" method
  idf <- log(1 + ncol(bmat) / Matrix::rowSums(bmat))
  # TF-IDF
  fast_tfidf = function(tf, idf) {
    tf <- t(tf)
    tf@x <- tf@x * rep.int(idf, diff(tf@p))
    tf <- t(tf)
    return(tf)
  }
  message("[info] fast log tf-idf")
  tf_idf_counts <- fast_tfidf(tf, idf)
  rownames(tf_idf_counts) <- rownames(bmat)
  colnames(tf_idf_counts) <- colnames(bmat)
  return(tf_idf_counts)
}
