#'
#' @title capOutlierQuantile
#' @description cap/normalize the trs with quantile
#' @param x a numeric vector
#' @param q_ceiling a numeric value indicating which quantile is used for ceiling the vector
#'
#' @return numeric vector
#' @export
#'
#' @importFrom stats quantile
#'
#' @examples
#' (x_ceiling <- capOutlierQuantile(x=1:10, q_ceiling=0.8))
#'
capOutlierQuantile <- function(x, q_ceiling=0.99){
  x[x>quantile(x, q_ceiling)] <- quantile(x, q_ceiling)
  return(x)
}
