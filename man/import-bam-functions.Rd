\name{import-bam-functions}
\alias{impBamGAL}
\alias{impBamGRL}


\title{Import Bam}
\description{
    Imports a specified genomic region from a bam file using a GenomicRanges object as search area. 
}
\usage{
  impBamGAL(UserDir, searchArea, XStag = FALSE, verbose = TRUE)
  impBamGRL(UserDir, searchArea, verbose = TRUE)
  }

\arguments{
  \item{UserDir}{The relative or full path of folder containing bam files.}
  \item{searchArea}{A \code{GenomicRanges object} that contains the regions of interest}
  \item{XStag}{Setting \code{XStag=TRUE} stores the strand specific information in the mcols slot 'XS'}
  \item{verbose}{Setting \code{verbose=TRUE} gives details of procedure during function run.}
}

\details{
These functions are wrappers to import bam files into R and store them into either GRanges, GAlignments or GappedAlignmentpairs objects. 

It is recommended to use the impBamGAL() which takes information of gaps into account. It is also possible to use the other variants as well, but then pre-filtering becomes important because gapped, intron-spanning reads will cause problems. This is because the GRanges objects can not handle if gaps are present and will then give a wrong result when calculating the allele (SNP) count table. 

If the sequence data is strand-specific you may want to set XStag=TRUE. The strand specific information will then be stored in the meta columns with column name 'XS'.

}

\value{
    \code{impBamGRL} returns a GRangesList object containing the RNA-seq reads in the region defined by the \code{searchArea} argument.
    \code{impBamGAL} returns a list with GAlignments objects containing the RNA-seq reads in the region defined by the \code{searchArea} argument.
    \code{funImpBamGAPL} returns a list with GappedAlignmentPairs object containing the RNA-seq reads in the region defined by the \code{searchArea} argument.
}

\note{
    A typical next step after the import of bam files is to obtain SNP information. This can be done either with the \code{\link{impBcfGRL}} and \code{\link{getAlleleCount}} functions. Alternatively the \code{\link{scanForHeterozygotes}} function provides R-based functionality for identifying heterozygote coding SNPs.
}

\seealso{
  \itemize{
    \item  The \code{\link{impBcfGRL}} for importing Bcf files.
  }
}

\author{Jesper R. Gadin, Lasse Folkersen}

\examples{
#Declare searchArea
searchArea <- GRanges(seqnames=c("17"), ranges=IRanges(79478301,79478361))

#Relative or full path  
pathToFiles <- system.file("extdata/ERP000101_subset", package="AllelicImbalance")

reads <- impBamGAL(pathToFiles,searchArea,verbose=FALSE)

}
\keyword{ bam }
\keyword{ import }


