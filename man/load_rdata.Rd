% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_rdata.R
\name{load_rdata}
\alias{load_rdata}
\title{\code{load_rdata}}
\usage{
load_rdata(fileName)
}
\arguments{
\item{fileName}{Name of the file to load.}
}
\value{
Data object.
}
\description{
Load processed data (\emph{.rda} format) using a function that assigns it
to a specific variable
(so you don't have to guess what the loaded variable name is).
}
\examples{
tmp <- tempfile()
save(mtcars, file = tmp)
mtcars2 <- load_rdata(tmp)
}
