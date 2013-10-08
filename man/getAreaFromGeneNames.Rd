\name{getAreaFromGeneNames}
\alias{getAreaFromGeneNames}
\title{Get Gene Area}
\description{
    Given a character vector with genesymbols and an OrgDb object, this function returns a GRanges giving the coordinates of the genes.
}
\usage{
	getAreaFromGeneNames(genesymbols, OrgDb, leftFlank=1000,
	rightFlank=1000, verbose=TRUE)}

\arguments{
  \item{genesymbols}{A character vector that contains genesymbols of genes from which we wish to retrieve the coordinates}
  \item{OrgDb}{An \code{OrgDb} object containing gene annotation}
  \item{leftFlank}{A \code{integer} specifying number of additional nucleotides before the genes}
  \item{rightFlank}{A \code{integer} specifying number of additional nucleotides after the genes}
  \item{verbose}{Setting \code{verbose=TRUE} makes function more talkative}
  
}
\details{
	This function is a convenience function that can be used to determine which genomic coordinates to specify to e.g. \code{impBamGAL} when retrieving reads. 
}
\value{
    \code{getAreaFromGeneNames} returns a GRanges object with genomic coordinates around the specified genes  
}

\author{Jesper R. Gadin, Lasse Folkersen}
\examples{
	#load example data
	data(ASEset)
	
	#get counts at the three positions specified in GRvariants
	library(org.Hs.eg.db )
	searchArea<-getAreaFromGeneNames(c("PAX8","TLR7"), org.Hs.eg.db)	
}

\keyword{ genes }
\keyword{ locations }


