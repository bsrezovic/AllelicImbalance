\name{getSnpIdFromLocation}
\alias{getSnpIdFromLocation}
\title{Get rsIDs from locations of SNP}
\description{
    Given a GRanges object of SNPs and a SNPlocs annotation, this function attempts to replace the names of the GRanges object entries with rs-IDs. 
}
\usage{
	getSnpIdFromLocation(GR, SNPloc, return.vector=FALSE, verbose=TRUE)}

\arguments{
  \item{GR}{A \code{GRanges} that contains positions of SNPs to look up}
  \item{SNPloc}{A \code{SNPlocs object} containing information on SNP locations (e.g. SNPlocs.Hsapiens.dbSNP.xxxxxxxx)}
  \item{return.vector}{Setting \code{return.vector=TRUE} returns vector with rsIds}
  \item{verbose}{Setting \code{verbose=TRUE} makes function more talkative}
}
\details{
	This function is used to try to identify the rs-IDs of SNPs in a GRanges object. 
}
\value{
    \code{getSnpIdFromLocation} returns the same GRanges object it was given with, but with updated with rs.id information. 
}

\author{Jesper R. Gadin, Lasse Folkersen}
\examples{
	#load example data
	data(ASEset)
	
	#get counts at the three positions specified in GRvariants
	library(SNPlocs.Hsapiens.dbSNP.20120608)
	updatedGRanges<-getSnpIdFromLocation(rowData(ASEset), SNPlocs.Hsapiens.dbSNP.20120608)
	rowData(ASEset)<-updatedGRanges

}

\keyword{ SNP }
\keyword{ rs-id }


