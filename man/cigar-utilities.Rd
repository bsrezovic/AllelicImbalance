\name{cigar-utilities}
\alias{cigar-utilities}
\alias{realCigarPosition}
\alias{realCigarPositions}
\alias{realCigarPositionsList}


\title{realCigarPosition}
\description{
   	From a GAlignments calculate the real corresponding position for each read based on its cigar.  
}
\usage{
    realCigarPosition(RleCigar, BpPos)
    realCigarPositions(RleCigar)
    realCigarPositionsList(RleCigarList)
}

\arguments{
  \item{RleCigar}{An \code{Rle} containing cigar information}
  \item{RleCigarList}{An \code{RleList} containing cigar information}
  \item{BpPos}{the absolute position on the chromosome of interest}
}
\details{
	The main intention for these functions are to be the internal functions for \code{scanForHeterozygotes} and \code{getAlleleCount}.
}
\value{
    \code{realCigarPosition} returns the new position
    \code{realCigarPositions} returns a vector with the corrected positions to be subsetted from a read.
    \code{realCigarPositionsList} returns a list where each element i a vector with the corrected positions to be subsetted from a read.
}

\seealso{
  \itemize{
    \item The \code{\link{scanForHeterozygotes}} which is a function to find possible heterozygote sites in a \code{\link[GenomicRanges]{GAlignmentsList}} object
  }
}

\author{Jesper R. Gadin}
\examples{
  RleCigarList <-  cigarToRleList("3M4I93M")
  BpPos <- 5

  newPos <- realCigarPosition(RleCigar=RleCigarList[[1]], BpPos)
  newPositions <- realCigarPositions(RleCigar=RleCigarList[[1]])
  newPositionsList <- realCigarPositionsList(RleCigarList=RleCigarList)
}

\keyword{ internal }


