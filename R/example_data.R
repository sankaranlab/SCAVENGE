#' Example data
#'
#' Example data of multiple data types.
#' @param name Name of the data to return.
#' @param return_path Whether to return just the path to the file (\code{TRUE})
#'  or the imported  R object itself \code{FALSE}.
#' @inheritParams gchromVAR::importBedScore
#' @export
#' @importFrom gchromVAR importBedScore
#' @importFrom SummarizedExperiment rowRanges
#' @examples
#' pbmc5krda <- example_data(name="pbmc5k_SE.rda")
example_data <- function(name=c("pbmc5k_SE.rda",
                                "mono.PP001.bed"),
                         return_path=FALSE,
                         colidx=5){
  name <- name[1]
  if(grepl("^pbmc5k_SE",name, ignore.case = TRUE)){
    path <- system.file('rda','pbmc5k_SE.rda', package='SCAVENGE')
    if(isTRUE(return_path)){
      return(path)
    } else {
      SE_pbmc5k <- load_rdata(fileName = path)
      return(SE_pbmc5k)
    }
  } else if (grepl("^mono.PP001",name, ignore.case = TRUE)){
    trait_file <- system.file('extdata','mono.PP001.bed', package='SCAVENGE')
    if(isTRUE(return_path)){
      return(path)
    } else {
      ### Recursion ####
      pbmc5krda <- example_data(name = "pbmc5k_SE.rda")
      trait_import <- gchromVAR::importBedScore(
        ranges = SummarizedExperiment::rowRanges(pbmc5krda),
        files = trait_file,
        colidx = colidx)
      return(trait_import)
    }
  #### Throw error ####
  } else {
    stp <- paste(
      "name must be one of:",
      paste("\n -",shQuote(eval(formals(example_data)$name)),collapse = "")
    )
    stop(stp)
  }
}
