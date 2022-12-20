# SCAVENGE 1.0.2

## New features

* Added a `NEWS.md` file to track changes to the package.
* Removed *docs*/*doc* folder in favor of *gh-pages* branch method.
* Renamed vignette "SCAVENGE-vignette.Rmd" --> "SCAVENGE.Rmd" 
  to make it the Get Started page for `pkgdown` site.
* Create *README.Rmd*
  - Add `rworkflows::use_badges()`
* Add *CITATION* file.
* New internal functions:
  - `load_rdata`
* New exported functions:
  - `example_data`
* Add `example_results` as *data.R*
* Add runnable example for all exported functions.
* Add *docker* vignette.

## Bug fixes

* Add *Remotes* to *DESCRIPTION*:
  - `gchromVAR`
  - `BuenColors`
* Add *URL* to *DESCRIPTION*.
* Remove duplicate `rmarkdown` entry in *DESCRIPTION*.
* Remove `lazyData: true` from *DESCRIPTION*.
* Update *.Rbuildignore*.
* Fix all: "no visible binding for global variable" CRAN Notes.
* Make necessary packages *IMPORT*s:
  - `SummarizedExperiment`
* Make argument explicit in vignettes.
* Replace `F` --> `FALSE`
* Make default `mycores=1`
* Replace all `1:n` with `seq_len(n)`
* `get_sigcell_simple`:
  - Create reproducible example
  - Remove defaults of variables that don't exist.
  - Remove "." usage and calling functions without the () syntax.
* Save output files as temp files by default.
* Improve vignettes
  - Use `scale_color_viridis_c` to simplify.
  - Reformat code to be more readable.
  - Add session info at the end.
* Remove `dplyr` as Dep in favor of native pipe.
  - Include as *Suggest* for vignette only



