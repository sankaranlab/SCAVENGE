#' @title cal_scalefactor
#'
#' @description Used for generate a scale factor from bias-corrected Z score
#' for TRS scale and normalization.
#'
#' @param z_score a numeric vector that is original Z scores for
#' individual cells calculated from gchromVAR.
#' @param percent_cut a numeric number (0, 1) indicating what percentage of
#' most enriched cells will be selected to calculate scale factor.
#'
#' @returns a numeric value
#' @export
#'
#' @examples
#' z_score <- example_results$mono_mat$z_score
#' sf <- cal_scalefactor(z_score = z_score)
cal_scalefactor <- function(z_score,
                            percent_cut=0.01){

  if(percent_cut<=0 | percent_cut>=1) {
    stop("percent_cut is a value between 0 and 1")
  }
  message("Scale factor is calculating from most enriched ",
          (100*percent_cut), "% of cells.")
  sf_idx <- rank(-z_score) <= floor(percent_cut*length(z_score))
  sf <- mean(z_score[sf_idx])
  return(sf)
}
