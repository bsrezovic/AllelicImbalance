\name{getDefaultMapBiasExpMean}
\alias{getDefaultMapBiasExpMean}
\alias{getDefaultMapBiasExpMean3D}

\title{Map Bias}
\description{
     an allele frequency list
}
\usage{
	getDefaultMapBiasExpMean(alleleCountList)
	getDefaultMapBiasExpMean3D(alleleCountList)
}
\arguments{
  \item{alleleCountList}{A \code{GRangesList object} containing read information}
}
\details{
	This function will assume there is no bias that comes from the mapping of reads, and therefore create a matrix with expected frequency of 0.5 for each allele.  

}
\value{
    \code{getDefaultMapBiasExpMean} returns a matrix with a default expected mean of 0.5 for every element.
}

\author{Jesper R. Gadin, Lasse Folkersen}
\examples{
	#load example data
	data(ASEset)
	#access SnpAfList
	alleleCountList <- alleleCounts(ASEset)
	#get default map bias exp mean
	matExpMean <- getDefaultMapBiasExpMean(alleleCountList)

}

\keyword{mapping bias}


