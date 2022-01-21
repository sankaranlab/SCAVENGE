#' @title get_sigcell_simple
#'
#' @description a function to determine statistically trait-enriched cell by permutation test
#'
#' @param knn_sparse_mat a sparse matrix used for network propagation, which indicates the adjacent matrix (m x m, where m is the cell number) of cell-to-cell network (M-kNN graph)
#' @param seed_idx a logical vector indicating seed cells (TRUE) and non-seed cells (FALSE) with length of m, where m is the cell number. The length and position are corresponding to knn_sparse_mat
#' @param topseed_npscore a numeric vector of real network propagation score
#' @param permutation_times an integer describe times of permutation test for each cell
#' @param true_cell_significance a numeric value between 0-1 indicating the significant threshold used to determine statistically trait-enriched cell
#' @param rda_output if output details of each permutation as a rda
#' @param out_rda if rda_output=T, an rda will be outputed with path and name specified
#' @param mycores how many cores used for permutation test
#' @param rw_gamma gamma from randomWalk_sparse function. Keep the same with what you used in the real setting
#'
#' @return a dataframe with two columns, the first one is seed index (the same with input), the second one is significance from permutation test. In addition, an R data is generated, which contains this dataframe and another dataframe depicting the network propagation score for each cell at each permuation
#' @export
#'
#' @import Matrix
#' @import parallel
#' @importFrom dplyr %>%
#'
#' @examples
#' \dontrun{
#' get_sigcell_simple(knn_sparse_mat=mutualknn30,
#' seed_idx=seed_p0.05,
#' permutation_times=1000,
#' true_cell_significance=0.05,
#' rda_output=F,
#' out_rda="true_cell_df.rda",
#'  mycores=4, rw_gamma=0.05)}
#'
get_sigcell_simple <- function(knn_sparse_mat=mutualknn30,
                               seed_idx=seed_p0.05,
                               topseed_npscore=topseed_npscore,
                               permutation_times=1000,
                               true_cell_significance=0.05,
                               rda_output=F,
                               out_rda="true_cell_df.rda",
                               mycores=4,
                               rw_gamma=0.05){
  if(permutation_times<100){
    warning("Permutation times less than 100")
  }
  stopifnot("true_cell_significance must be a numeric value between 0, 1" = (true_cell_significance>0 & true_cell_significance<1))
  message("Get started!")
  cell_mat <- data.frame(cell=1:nrow(knn_sparse_mat), degree=colSums(knn_sparse_mat))
  cell_table <- data.frame(table(cell_mat$degree))

  seed_mat_top  <- data.frame(seed=which(seed_idx), degree=colSums(knn_sparse_mat[, seed_idx]))
  # summary(seed_mat_top $degree)
  seed_table_top <- data.frame(table(seed_mat_top$degree))
  xx_top <- tapply(cell_mat[, 1], cell_mat[, 2], list)
  xx2_top <- xx_top[names(xx_top) %in% seed_table_top$Var1]
  # permutation_score_top <- data.frame(matrix(nrow=nrow(knn_sparse_mat), ncol=1000))
    permutation_score_top <- mclapply(1:permutation_times, mc.cores = mycores, function(i){
      sampled_cellid <- xx2_top %>%
        mapply(sample, ., seed_table_top$Freq) %>%
        unlist %>%
        sort
      xx <- randomWalk_sparse(intM=knn_sparse_mat, queryCells=rownames(knn_sparse_mat)[as.numeric(sampled_cellid)], gamma=rw_gamma)
      if (i %% 100 == 0) {message(i)}
      return(xx)
    }
    )

  names(permutation_score_top) <- paste0("permutation_", 1:permutation_times)
  permutation_score_top <- data.frame(sapply(permutation_score_top, c))
  # permutation_score_top <- do.call(cbind.data.frame, permutation_score_top)
  permutation_df_top <- data.frame(matrix(nrow=nrow(knn_sparse_mat), ncol=permutation_times))

  permutation_df_top <- apply(permutation_score_top, 2, function(x) { temp <- x > topseed_npscore; return(temp) } )
  message("cells passed 0.001 threshold: ", round(sum(rowSums(permutation_df_top) <= 0.001*permutation_times)*100/nrow(permutation_df_top), 2), "%")
  message("cells passed 0.01 threshold: ", round(sum(rowSums(permutation_df_top) <= 0.01*permutation_times)*100/nrow(permutation_df_top), 2), "%")
  message("cells passed 0.05 threshold: ", round(sum(rowSums(permutation_df_top) <= 0.05*permutation_times)*100/nrow(permutation_df_top), 2), "%")
  true_cell_top_idx <- rowSums(permutation_df_top) <= true_cell_significance*permutation_times
  message("your emprical P value threshold: ", true_cell_significance)
  message("what propertion of enriched cells over all cells: ", round((sum(true_cell_top_idx)*100)/nrow(permutation_df_top), 2), "%") # how many propertion true cell over all cells
  message("what propertion of seed cells that are enriched cells: ",  round(sum(true_cell_top_idx & seed_idx)*100/sum(seed_idx), 2), "%") # how many propertion of seed were true cells
  message("fold of true cell over seed: ", round(sum(true_cell_top_idx)/sum(seed_idx), 2)) # fold of true cell over seed

  true_cell_top_filter_idx <-
    ture_cell_df <- data.frame(seed_idx,
                               true_cell_top_idx)
  if(rda_output){
    save(ture_cell_df, permutation_score_top, file=out_rda)
  }
  return(ture_cell_df)
}
