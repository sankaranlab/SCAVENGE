#' \code{load_rdata}
#'
#' Load processed data (\emph{.rda} format) using a function that assigns it
#' to a specific variable
#'  (so you don't have to guess what the loaded variable name is).
#'
#' @param fileName Name of the file to load.
#'
#' @return Data object.
#'
#' @export
#' @examples
#' tmp <- tempfile()
#' save(mtcars, file = tmp)
#' mtcars2 <- load_rdata(tmp)
load_rdata <- function(fileName) {
    load(fileName)
    get(ls()[ls() != "fileName"])
}
