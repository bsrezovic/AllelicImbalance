\name{coverageMatrixListFromGAL}
\alias{coverageMatrixListFromGAL}


\title{coverage matrix of GAlignmentsList}
\description{
  Get coverage per nucleotide for reads covering a region
}

\usage{
  coverageMatrixListFromGAL(BamList,strand=NULL,ignore.empty.bam.row=TRUE)
}

\arguments{
  \item{BamList}{GAlignmentsList containing reads over the region to calculate coverage }
  \item{strand}{strand has to be '+' or '-'}
  \item{ignore.empty.bam.row}{argument not in use atm }
}

\details{
  a convenience function to get the coverage from a list of reads stored in GAlignmnetsList, and returns by default a list with one matrix, and information about the genomic start and stop positions.
}

\author{Jesper R. Gadin}

\examples{
	r <- reads
	seqlevels(r) <- "17"
	covMatList <- coverageMatrixListFromGAL(BamList=r, strand="+")
}
\keyword{ coverage }


