#' @title max_min_scale
#'
#' @description a function for max-min normalization for
#' network propagation score.
#'
#' @param x a numeric vector
#'
#' @returns a numeric vector
#' @export
#'
#' @examples
#' x <- rnorm(10)
#' scale_vec <- max_min_scale(x = x)
max_min_scale <- function(x){
  scale_vec <- (x-min(x))/(max(x)-min(x))
  return(scale_vec)
}
