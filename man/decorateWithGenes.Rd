\name{decorateWithGenes}
\alias{decorateWithGenes}

\title{decorateWithGenes}
\description{
   	Internal function that can draw gene regions on pre-specified surfaces. Necessary for the genomic-location plots.  
}
\usage{
	decorateWithGenes(x,genesInRegion,xlim,ylim,chromosome)
}

\arguments{
  \item{x}{\code{ASEset} object}
  \item{genesInRegion}{\code{GRanges} object with gene regions. Can be obtained using \code{getGenesFromAnnotation}}
  \item{xlim}{xlim values for the pre-specified surface}
  \item{ylim}{ylim values for the pre-specified surface}
  \item{chromosome}{character}
}
\details{
	The main intention of this function is to be used when plotting several bar plots in the same window. This function add gene regions under the bars. 
}
\value{
    \code{decorateWithGenes} returns nothing, but draws genes
}
\seealso{
  \itemize{
    \item The \code{\link{locationplot}} which is uses this function internally.
    \item The \code{\link{decorateWithExons}} which is another similar function that \code{\link{locationplot}} uses internally. 
  }
}



\author{Jesper R. Gadin, Lasse Folkersen}
\examples{
  data(ASEset)
}

\keyword{ internal }


