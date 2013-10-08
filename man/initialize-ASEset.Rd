
\name{initialize-ASEset}
\alias{initialize-ASEset}
\alias{ASEsetFromCountList}


\title{Initialize ASEset}
\description{ Functions to construct ASEset objects }
\usage{
ASEsetFromCountList(rowData, countListNonStranded = NULL, countListPlus =
                 NULL, countListMinus = NULL, countListUnknown = NULL,
                 colData = NULL, mapBiasExpMean = NULL, verbose =
                 FALSE, ...)	
}

\arguments{
  \item{rowData}{A \code{GenomicRanges object} that contains the variants of interest}
  \item{countListNonStranded}{A \code{list} where each entry is a matrix with allele counts as columns and sample counts as rows}
  \item{countListPlus}{A \code{list} where each entry is a matrix with allele counts as columns and sample counts as rows}
  \item{countListMinus}{A \code{list} where each entry is a matrix with allele counts as columns and sample counts as rows}
  \item{countListUnknown}{A \code{list} where each entry is a matrix with allele counts as columns and sample counts as rows}
  \item{colData}{A \code{DataFrame} object containing sample specific data}
  \item{mapBiasExpMean}{A 3D \code{array} where the SNPs are in the 1st dimension, samples in the 2nd dimension and variants in the 3rd dimension.}
  \item{verbose}{Makes function more talkative}
  \item{...}{arguments passed on to SummarizedExperiment constructor}
}

\details{
	The resulting \code{ASEset} object is based on the \code{SummarizedExperiment}, and will therefore inherit the same accessors and ranges operations.

	countListNonStranded, countListPlus, countListMinus and countListUnknown are i.e. the outputs from the \code{\link{getAlleleCount}} function. 
 }

 \note{
	\code{ASEsetFromCountList} requires the same input data as an SummarizedExperiment, but with minimum one assay for the allele counts.
 }

\value{
    \code{ASEsetFromCountList} returns an \code{ASEset} object.
}

\author{Jesper R. Gadin, Lasse Folkersen}

\seealso{
	\itemize{
		\item The \code{\link[GenomicRanges]{SummarizedExperiment}} for ranges operations.
	}
}

\examples{
#make example alleleCountListPlus
set.seed(42)
countListPlus <- list()
snps <- c("snp1","snp2","snp3","snp4","snp5")
for(snp in snps){
	count<-matrix(rep(0,20),ncol=5,dimnames=list(
		c("sample1","sample2","sample3","sample4"),
		c("A","T","G","C","del")))
	#insert random counts in two of the alleles 
	for(allele in sample(c("A","T","G","C","del"),2)){
		count[,allele]<-as.integer(rnorm(4,mean=50,sd=10))
	}
	countListPlus[[snp]] <- count
}

#make example alleleCountListMinus
countListMinus <- list()
snps <- c("snp1","snp2","snp3","snp4","snp5")
for(snp in snps){
	count<-matrix(rep(0,20),ncol=5,dimnames=list(
		c("sample1","sample2","sample3","sample4"),
		c("A","T","G","C","del")))
	#insert random counts in two of the alleles 
	for(allele in sample(c("A","T","G","C","del"),2)){
		count[,allele]<-as.integer(rnorm(4,mean=50,sd=10))
	}
	countListMinus[[snp]] <- count
}


#make example rowData
rowData <- GRanges(
		seqnames = Rle(c("chr1", "chr2", "chr1", "chr3", "chr1")),
         ranges = IRanges(1:5, width = 1, names = head(letters,5)),
         snp = paste("snp",1:5,sep="")
         )
#make example colData
colData <- DataFrame(Treatment=c("ChIP", "Input","Input","ChIP"), 
 	row.names=c("ind1","ind2","ind3","ind4"))

#make ASEset 
a <- ASEsetFromCountList(rowData, countListPlus=countListPlus, 
countListMinus=countListMinus, colData=colData)

}
\keyword{ ASEset }
\keyword{ ASEsetFromCountList }
