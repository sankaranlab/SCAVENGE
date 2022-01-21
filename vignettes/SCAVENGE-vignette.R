## ---- message = FALSE, warning = FALSE----------------------------------------
library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(igraph)

set.seed(9527)

## ----message=TRUE, warning=FALSE, cache=FALSE---------------------------------
trait_file <- paste0(system.file('extdata', package='SCAVENGE'), "/mono.PP001.bed")
pbmc5krda <- paste0(system.file('rda', package='SCAVENGE'), "/pbmc5k_SE.rda")
load(pbmc5krda)

