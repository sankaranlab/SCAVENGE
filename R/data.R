#' Example results
#'
#' @description
#' Example results for SCAVENGE.
#' \code{
#' trait_import <- example_data(name="mono.PP001.bed")
#' SE_pbmc5k <- example_data(name="pbmc5k_SE.rda")
#' genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#' SE_pbmc5k <- addGCBias(object = SE_pbmc5k,
#'                        genome = genome)
#' SE_pbmc5k_bg <- getBackgroundPeaks(object = SE_pbmc5k,
#'                                    niterations = 200)
#' SE_pbmc5k_DEV <- computeWeightedDeviations(object = SE_pbmc5k,
#'                                            weights = trait_import,
#'                                            background_peaks = SE_pbmc5k_bg)
#' z_score_mat <- data.frame(colData(SE_pbmc5k),
#'                           z_score=c(t(assays(SE_pbmc5k_DEV)[["z"]])) )
#' scale_factor <- cal_scalefactor(z_score = z_score_mat$z_score,
#'                                 percent_cut = 0.01)
#' seed_idx <- seedindex(z_score_mat$z_score, 0.05)
#' peak_by_cell_mat <- SummarizedExperiment::assay(SE_pbmc5k)
#' tfidf_mat <- tfidf(bmat=peak_by_cell_mat)
#' lsi_mat <- do_lsi(mat = tfidf_mat, dims = 30)
#' mutualknn30 <- getmutualknn(lsimat = lsi_mat, num_k = 30)
#' np_score <- randomWalk_sparse(intM=mutualknn30,
#'                               queryCells = rownames(mutualknn30)[seed_idx],
#'                               gamma=0.05)
#' omit_idx <- np_score==0
#' mutualknn30 <- mutualknn30[!omit_idx, !omit_idx]
#' np_score <- np_score[!omit_idx]
#' TRS <- capOutlierQuantile(np_score, 0.95) |> max_min_scale()
#' TRS <- TRS * scale_factor
#' mono_mat <- data.frame(z_score_mat[!omit_idx, ],
#'                        seed_idx[!omit_idx],
#'                        np_score,
#'                        TRS)
#' #### Merge into one list ####
#' example_results <- list(seed_idx=seed_idx,
#'                         # tfidf_mat=tfidf_mat,
#'                         # lsi_mat=lsi_mat,
#'                         mutualknn30=mutualknn30,
#'                         mono_mat=mono_mat)
#' usethis::use_data(example_results, overwrite = TRUE)
#' }
#' @format List
#' @usage data("example_results")
"example_results"
