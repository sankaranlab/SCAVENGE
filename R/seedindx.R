#' @title seedindx
#'
#' @description used for generate seed cell index from bias-corrected Z score
#'
#' @param z_score a numeric vector that is original Z scores for individual
#' cells calculated from gchromVAR.
#' @param percent_cut a numeric number (0, 1) indicating what percentage
#' of cells will be selected as seed cells if too many cells fit
#' the P value<0.05 cutoff.
#'
#' @returns a logical vector indicating seed cells or not.
#' @export
#'
#' @importFrom stats pnorm
#'
#' @examples
#' z_score <- example_results$mono_mat$z_score
#' seed_idx <- seedindex(z_score = z_score, percent_cut = 0.1)
seedindex <- function(z_score,
                      percent_cut=0.05){

  if(percent_cut<=0 | percent_cut>=1) {
    stop("percent_cut is a value between 0 and 1")
  }
  z2pvalue <- stats::pnorm(q = z_score,
                           lower.tail = FALSE)
  message("Cells with enriched P < 0.05: ", sum(z2pvalue<=0.05))
  s_percent <- round(sum(z2pvalue<=0.05)*100/length(z2pvalue), 2)
  message("Percent: ", s_percent, "%")
  if(s_percent>(100*percent_cut)){
    message("The top ", (100*percent_cut),
            "% of cells ", "(N=", floor(percent_cut*length(z2pvalue)),
            ") were selected as seed cells")
    seed_idx <- rank(-z_score) <= floor(percent_cut*length(z2pvalue))
  } else {
    message("These cells were selected as seed cells")
    seed_idx <- z2pvalue<=0.05
  }
  return(seed_idx)
}
