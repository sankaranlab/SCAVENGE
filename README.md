# SCAVENGE

### Overview:

Co-localization approaches using genetic variants and single-cell epigenomic data are unfortunately uninformative for many cells given the extensive sparsity across single-cell profiles. Therefore, only a few cells from the truly relevant population demonstrate reliable phenotypic relevance. The global high-dimensional features of individual single cells are sufficient to represent the underlying cell identities or states, which enables the relationships among such cells to be readily inferred. By taking advantage of these attributes, SCAVENGE identifies the most phenotypically-enriched cells by co-localization and explores the transitive associations across the cell-to-cell network to assign each cell a probability representing the cell’s relevance to those phenotype-enriched cells via network propagation.

We developed a novel enrichment method (**SCAVENGE**) (Single Cell Analysis of Variant Enrichment through Network propagation of GEnomic data) that can discriminate between closely related cell types/states and score single cells for GWAS enrichment. 


<div align=center> <img src="image/schematic-view_1.png" width="680" height="278"> </div> 

<p align="center">Schematic view of SCAVENGE</p>  



We've implemented **SCAVENGE** as an `R` package for computing single-cell based GWAS enrichments from fine-mapped posterior probabilities and quantitative epigenomic data (i.e. scATAC-seq and potentially other single-cell epigenome profiling methods).  
As single-cell genomic datasets grow in volume, we expect SCAVENGE will have great promise for efficiently uncovering relevant cell populations for more phenotypes or functions in different scenarios, which may expand beyond the complex trait genetic variants we have examined here. We welcome you to use SCAVENGE to discover more phenotype relevant cells!

### Installation:

Once all of the dependencies for `SCAVENGE` are installed, the package can be installed 
directly from GitHub by typing the following into an `R` console:

```
devtools::install_github("https://github.com/sankaranlab/SCAVENGE")
```
### Tutorial:
This web resource and vignette compiliation shows how to reproduce results of SCAVENGE analysis with monocyte count on a 10X PBMC dataset [[Vignette-pdf]](doc/SCAVENGE-vignette.pdf), [[Vignette-R markdown code]](doc/SCAVENGE-vignette.Rmd). 


### Citation:
If you used or adapted SCAVENGE in your study, please cite our paper [[*Nat Biotechnol*]](https://www.nature.com/articles/s41587-022-01341-y)  || [[*PubMed*]](https://pubmed.ncbi.nlm.nih.gov/35668323/).   
*Variant to function mapping at single-cell resolution through network propagation.*


### Contact:
If you run into issues and would like to report them, you can use the "Issues" tab on the left hand side.  
Alternatively, you can contact authors: fyu{at}broadinstitute.org or sankaran{at}broadinstitute.org.  

