<a href="http://www.bioconductor.org/packages/devel/bioc/html/AllelicImbalance.html#since"><img border="0" src="http://www.bioconductor.org/shields/years-in-bioc/AllelicImbalance.svg" title="How long since the package was first in a develd Bioconductor version (or is it in devel only)."></a> <a href="http://bioconductor.org/packages/stats/bioc/AllelicImbalance.html"><img border="0" src="http://www.bioconductor.org/shields/downloads/AllelicImbalance.svg" title="Percentile (top 5/20/50% or 'available') of downloads over last 6 full months. Comparison is done across all package categories (software, annotation, experiment)."></a> <a href="https://support.bioconductor.org/t/AllelicImbalance/"><img border="0" src="http://www.bioconductor.org/shields/posts/AllelicImbalance.svg" title="Support site activity, last 6 months: tagged questions/avg. answers per question/avg. comments per question/accepted answers, or 0 if no tagged posts."></a> <a href="http://www.bioconductor.org/packages/devel/bioc/html/AllelicImbalance.html#svn_source"><img border="0" src="http://www.bioconductor.org/shields/commits/bioc/AllelicImbalance.svg" title="average Subversion commits (to the devel branch) per month for the last 6 months"></a>

Status: Travis CI [![Build Status](https://travis-ci.org/lcolladotor/AllelicImbalance.svg?branch=master)](https://travis-ci.org/lcolladotor/AllelicImbalance)

Bioc-release <a href="http://www.bioconductor.org/packages/release/bioc/html/AllelicImbalance.html#archives"><img border="0" src="http://www.bioconductor.org/shields/availability/release/AllelicImbalance.svg" title="Whether the package is available on all platforms; click for details."></a> <a href="http://bioconductor.org/checkResults/release/bioc-LATEST/AllelicImbalance/"><img border="0" src="http://www.bioconductor.org/shields/build/release/bioc/AllelicImbalance.svg" title="build results; click for full report"></a>

Bioc-devel <a href="http://www.bioconductor.org/packages/devel/bioc/html/AllelicImbalance.html#archives"><img border="0" src="http://www.bioconductor.org/shields/availability/devel/AllelicImbalance.svg" title="Whether the package is available on all platforms; click for details."></a> <a href="http://bioconductor.org/checkResults/devel/bioc-LATEST/AllelicImbalance/"><img border="0" src="http://www.bioconductor.org/shields/build/devel/bioc/AllelicImbalance.svg" title="build results; click for full report"></a>

The Bioc-release and Bioc-devel builds are updated every day at 10:30AM Seattle
time.

The latest bioconductor releases and the bioconductor development branch can be
found either on [github][bioc-mirror-branches] or at the home of
[bioconductor][bioc-home], 
* [release][bioc-release]
* [devel][bioc-devel]

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
install_github("AllelicImbalance","pappewaio")
```

Alternatively, you can easily install the release version of **AllelicImbalance** from
Bioconductor.

```r
source("http://bioconductor.org/biocLite.R")
biocLite("AllelicImbalance")
```

## Citation

Jesper R. Gådin, Ferdinand M. van’t Hooft, Per Eriksson and Lasse Folkersen.
(2015). AllelicImbalance: an R/ bioconductor package for detecting, managing,
and visualizing allele expression imbalance data from RNA sequencing. BMC
Bioinformatics. [link](http://www.biomedcentral.com/1471-2105/16/194)

## Bug Reports / Issues

The official development site for `AllelicImbalance` is at
https://github.com/pappewaio/AllelicImbalance. For bug-reports / issues /
comments regarding the development branch you may preferably submit there, and
for the release version it might be more appropriate to submit at the
bioconductor support-site which can be found [HERE][bioc-support-site]


[vignette-release]: http://bioconductor.org/packages/release/bioc/vignettes/AllelicImbalance/inst/doc/AllelicImbalance-vignette.pdf "AllelicImbalance Vignette"
[vignette-devel]: http://bioconductor.org/packages/devel/bioc/vignettes/AllelicImbalance/inst/doc/AllelicImbalance-vignette.pdf "AllelicImbalance Vignette"
[bioc-mirror-branches]: https://github.com/Bioconductor-mirror/AllelicImbalance/branches "bioc-mirror-branches"
[bioc-home]: www.bioconductor.org "Bioconductor home"
[bioc-release]: http://bioconductor.org/packages/release/bioc/html/AllelicImbalance.html "Bioconductor release"
[bioc-devel]: http://bioconductor.org/packages/devel/bioc/html/AllelicImbalance.html "Bioconductor devel"
[bioc-support-site]: https://support.bioconductor.org/t/AllelicImbalance "support site AllelicImbalance"
[devtoolsLink]: https://github.com/hadley/devtools "devtools"

