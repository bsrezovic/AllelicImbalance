\name{ASEset-class}
\docType{class}

%ASEset-class
\alias{ASEset-class}

%constructors
\alias{ASEset}

%\alias{ASEset,ASEset-method}

%accessors
\alias{alleleCounts}
\alias{mapBias}
\alias{fraction}

\alias{alleleCounts,ASEset-method}
\alias{mapBias,ASEset-method}
\alias{fraction,ASEset-method}


\title{ASEset objects}
\description{
Object that holds allele counts, genomic positions and map-bias for a set of SNPs
}
\section{Constructor}{
    ASEsetFromCountList(rowData, countListNonStranded = NULL, countListPlus = NULL,
                 countListMinus = NULL, countListUnknown = NULL,
                 colData = NULL, mapBiasExpMean = array(),verbose=FALSE
                 ...)	
	
	\describe{
		
		Arguments:
		\item{rowData}{A \code{GenomicRanges object} that contains the variants of interest}
  		\item{countListNonStranded}{A \code{list} where each entry is a matrix with allele counts as columns and sample counts as rows}
  		\item{countListPlus}{A \code{list} where each entry is a matrix with allele counts as columns and sample counts as rows}
		\item{countListMinus}{A \code{list} where each entry is a matrix with allele counts as columns and sample counts as rows}
		\item{countListUnknown}{A \code{list} where each entry is a matrix with allele counts as columns and sample counts as rows}
		\item{colData}{A \code{DataFrame} object containing sample specific data}
		\item{mapBiasExpMean}{A 3D \code{array} describing mapping bias. The SNPs are in the 1st dimension, samples in the 2nd dimension and variants in the 3rd dimension.}
		\item{verbose}{Makes function more talkative}
		\item{...}{arguments passed on to SummarizedExperiment constructor}
	}
}
\section{Accessors}{
  In the following code snippets, \code{x} is an ASEset object.

	\describe{
		\item{}{
		\code{alleleCounts(x, strand="nonStranded")}:
		Get SNP allele count list. Strand options are "nonStranded", "+","-" or "*".
		}
		\item{}{
		\code{mapBias(x)}:
		Get SNP mapping bias matrix. 
		}
		\item{}{
		\code{fraction(x, strand="nonStranded")}:
		Get SNP allele fraction matrix. strand options are "nonStranded", "+","-" or"*".
		}
	}
}


\section{Subsetting}{
  Subset like a \code{\link[GenomicRanges]{SummarizedExperiment}}.

}

\value{An object of class ASEset containing location information and allele counts for a number of SNPs measured in a number of samples on various strand, as well as mapBias information. All data is stored in a manner similar to the \code{\link[GenomicRanges]{SummarizedExperiment}} class.}

\details{
	An ASEset object differs from a regular SummarizedExperiment object in that the assays contains an array instead of matrix. This array has ranges on the rows, sampleNames on the columns and variants in the third dimension.

	It is possible to use the commands barplot and locationplot on an ASEset object see more details in \code{\link{barplot}} and \code{\link{locationplot}}.

	Four different alleleCount options are available. The simples one is the nonStranded option, and is experiments where the strand information is not known. In this option both the plus, minus and unknown strand will be counted and present. The unknown strand is when the aligner could not find any strand associated with the read. Then there are an option too add plus and minus stranded data. When using this, it is essential to make sure that the RNA-seq experiment under analysis has in fact been created so that correct strand information was obtained. 

} 

\author{Jesper R. Gadin, Lasse Folkersen}

\seealso{
	\itemize{
		\item The \code{\link[GenomicRanges]{SummarizedExperiment}} for ranges operations.
	}
}

\examples{

	#make example countList
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
	colData=colData)


}


\keyword{class}
\keyword{ASEset}

