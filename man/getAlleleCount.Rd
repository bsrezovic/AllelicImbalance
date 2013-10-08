\name{getAlleleCount}
\alias{getAlleleCount}
\title{snp count data}
\description{
    Given the positions of known SNPs, this function returns allele counts from a BamGRL object
}
\usage{
getAlleleCount(BamList, GRvariants, strand, verbose=TRUE)}

\arguments{
 \item{BamList}{A \code{GAlignmentsList object} or \code{GRangesList object} containing data imported from a bam file}
  \item{GRvariants}{A \code{GRanges object} that contains positions of SNPs to retrieve}
  \item{strand}{A length 1 \code{character} with value 'nonStranded', '+', '-', or '*'.  This argument determines if \code{getAlleleCount} will retrieve counts from all reads, or only from reads marked as '+', '-' or '*' (unknown), respectively.}
  \item{verbose}{Setting \code{verbose=TRUE} makes function more talkative}
}
\details{
	This function is used to retrieve the allele counts from specified positions in a set of RNA-seq reads. The \code{BamList} argument will typically have been created using the \code{impBamGAL} function on bam-files. The \code{GRvariants} is either a GRanges with user-specified locations or else it is generated through scanning the same bam-files as in \code{BamList} for heterozygote locations (e.g. using \code{scanForHeterozygotes}). The GRvariants will currently only accept locations having width=1, corresponding to SNPs.  In the \code{strand} argument, specifying 'nonStranded' is the same as retrieving the sum count of '+' and '-' reads (and '*' unknown in case these are found in the bam file). 'nonStranded' is the default behaviour and can be used when the RNA-seq experiments strand information is not available. 
}
\value{
    \code{getAlleleCount} returns a list of several data.frame objects, each storing the count data for one SNP.  
}

\author{Jesper R. Gadin, Lasse Folkersen}

\seealso{
  \itemize{
    \item The \code{\link{scanForHeterozygotes}} which is a function to find possible heterozygote sites in a \code{\link[GenomicRanges]{GAlignmentsList}} object
  }
}

\examples{
	#load example data
	data(reads)
	data(GRvariants)

	#get counts at the three positions specified in GRvariants
	alleleCount <- getAlleleCount(BamList=reads,GRvariants,
		strand="nonStranded")
	
	#if the reads had contained stranded data, these two calls would 
	#have given the correct input objects for getAlleleCount
	alleleCountPlus <- getAlleleCount(BamList=reads,GRvariants,
		strand="+")
	alleleCountMinus <- getAlleleCount(BamList=reads,GRvariants,
		strand="-")
	
}

\keyword{ SNP }
\keyword{ count }


