## ---- message = FALSE, warning = FALSE----------------------------------------
devtools::load_all()
library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
set.seed(9527)

## ----message=TRUE, warning=FALSE, cache=FALSE---------------------------------
trait_file <- paste0(system.file('inst/extdata', package='SCAVENGE'), "/mono.PP001.bed")
pbmc5krda <- paste0(system.file('inst/rda', package='SCAVENGE'), "/pbmc5k_SE.rda")
load(pbmc5krda)

