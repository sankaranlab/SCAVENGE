#' @title max_min_scale
#'
#' @description a function for max-min normalization for network propagation score
#'
#' @param x a numeric vector
#'
#' @return a numeric vector
#' @export
#'
#' @examples
#' \dontrun{
#' trs <- max_min_scale(ceiling_np)}
#'
max_min_scale <- function(x){
  scale_vec <- (x-min(x))/(max(x)-min(x))
  return(scale_vec)
}
