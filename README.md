
This is a fork of the AllelicImbalance package (original dev page at https://github.com/pappewaio/AllelicImbalance).
It solves some bugs that cropped up due to the package not being updated in a while. All of these are in functions that
 retrieve gene names / manipulate tabular data. Core allele imbalance analysis functions remain unchanged.
If used please cite the original authors cited below.

# AllelicImbalance

A Bioconductor package for allelic imbalance analysis on RNA sequencing data.
The aim is to support infrastructure for storage of summarized allele counts and
functions to quickly retrieve important figures for AI analyses. The package
also provides SNP or region specific visualization and interaction with other
common bioconductor functionality.

See the "DESCRIPTION" file for additional requirements.

## Documentation

The vignette provides an introduction to the allelic imbalance
concept and show examples of how the AllelicImbalance package can be used to
find AI-events and measure the global or local mapping bias from the alignment
of RNA-seq reads to a reference genome. 

* [vignette-release][vignette-release]
* [vignette-devel][vignette-devel]

## Installation

Installation of the devel version of this package from Github requires the
[devtools][devtoolsLink] package.

```r
install.packages("devtools")
library(devtools)
install_github("bsrezovic/AllelicImbalance")
```


## Citation

Jesper R. Gådin, Ferdinand M. van’t Hooft, Per Eriksson and Lasse Folkersen.
(2015). AllelicImbalance: an R/ bioconductor package for detecting, managing,
and visualizing allele expression imbalance data from RNA sequencing. BMC
Bioinformatics. [link](http://www.biomedcentral.com/1471-2105/16/194)



[vignette-release]: http://bioconductor.org/packages/release/bioc/vignettes/AllelicImbalance/inst/doc/AllelicImbalance-vignette.pdf "AllelicImbalance Vignette"
[vignette-devel]: http://bioconductor.org/packages/devel/bioc/vignettes/AllelicImbalance/inst/doc/AllelicImbalance-vignette.pdf "AllelicImbalance Vignette"
[bioc-mirror-branches]: https://github.com/Bioconductor-mirror/AllelicImbalance/branches "bioc-mirror-branches"
[bioc-home]: https://www.bioconductor.org/ "Bioconductor home"
[bioc-release]: http://bioconductor.org/packages/release/bioc/html/AllelicImbalance.html "Bioconductor release"
[bioc-devel]: http://bioconductor.org/packages/devel/bioc/html/AllelicImbalance.html "Bioconductor devel"
[bioc-support-site]: https://support.bioconductor.org/t/AllelicImbalance "support site AllelicImbalance"
[devtoolsLink]: https://github.com/hadley/devtools "devtools"

