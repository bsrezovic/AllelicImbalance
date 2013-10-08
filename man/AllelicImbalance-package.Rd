\name{AllelicImbalance-package}
\alias{AllelicImbalance-package}
\alias{AllelicImbalance}
\docType{package}
\title{A package meant to provide all basic functions for high-throughput allele specific expression analysis}

\description{
Package \code{AllelicImbalance} has functions for importing, filtering and plotting high-throughput data to make an allele specific expression analysis. A major aim of this package is to provide functions to collect as much information as possible from regions of choice, and to be able to explore the allelic expression of that region in detail.
}

\details{
	\tabular{ll}{
	Package: \tab AllelicImbalance\cr
	Type: \tab Package\cr
	Version: \tab 0.99.0\cr
	Date: \tab 2013-08-27\cr
	License: \tab GPL-3\cr
	}
}


\author{
Author: Jesper Robert Gadin
Author: Lasse Folkersen

Maintainer: Jesper Robert Gadin <j.r.gadin@gmail.com>

}
\references{
Reference to published application note (work in progress)

}
\section{Overview - standard procedure}{
Start out creating a \code{GRange} object defining the region of interest. This can also be done using \code{getAreaFromGeneNames} providing gene names as arguments. Then use \code{BamImpGAList} to import reads from that reagion and find potential SNPs using \code{scanForHeterozygotes}. Then retrieve the allele counts of heterozygote sites by the function \code{getAlleleCount}. With this data create an \code{ASEset}. At this point all pre-requisites for a 'basic' allele specific expression analysis is available. Two ways to go on could be to apply \code{\link{chisq.test}} or \code{\link{barplot}} on this ASEset object.
}


\seealso{
	\itemize{
		\item code{?ASEset}
	}
}

\keyword{ package }