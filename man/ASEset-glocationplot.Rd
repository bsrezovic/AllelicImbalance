\name{ASEset-locationplot}
\docType{methods}
\alias{ASEset-locationplot}
\alias{locationplot}
\alias{locationplot,ASEset-method}


\title{locationplot ASEset objects}
\description{
	plotting ASE effects over a specific genomic region
}
\usage{
\S4method{locationplot}{ASEset}(
	x,
	type="fraction",
	strand="nonStranded",
	yaxis=TRUE,
	xaxis=FALSE,
	xlab=FALSE,
	ylab=TRUE,
	legend.colnames = "", 
	size=c(0.8,1),
	main=NULL,
	pValue=FALSE,
	cex.main=0.7,
	cex.ylab=0.6,
	cex.legend= 0.5,
	OrgDb=NULL,
	TxDb=NULL,
	verbose=TRUE,
	...)

}
\arguments{
	\item{x}{ an ASEset object.}
	\item{type}{"fraction" or "count"}
	\item{strand}{"+","-","*","both" or "nonStranded". This argument determines which strand is plotted. See \code{getAlleleCounts} for more information on strand.}
	\item{yaxis}{wheter the y-axis is to be displayed or not}
	\item{xaxis}{wheter the x-axis is to be displayed or not}
	\item{ylab}{showing labels for the tic marks}
	\item{xlab}{showing labels for the tic marks}
	\item{legend.colnames}{gives colnames to the legend matrix}
	\item{size}{will give extra space in the margins of the inner plots}
	\item{main}{text to use as main label}
	\item{pValue}{Display p-value}
	\item{cex.main}{set main label size}
	\item{cex.ylab}{set ylab label size}
	\item{cex.legend}{set legend label size}
	\item{OrgDb}{an OrgDb object from which to plot a gene map. If given together with argument TxDb this will only be used to extract genesymbols.}
	\item{TxDb}{a TxDb object from which to plot an exon map.}
	\item{verbose}{Setting \code{verbose=TRUE} gives details of procedure during function run}
	\item{...}{arguments passed on to barplot function}
}	


\details{
	The locationplot methods visualises how fractions are distributed over a larger region of genes on one chromosome. It takes and ASEset object as well as additional information on plot type (see \code{\link{barplot}}), strand type (see \code{\link{getAlleleCounts}}), colouring, as well as annotation. The annotation is taken either from the bioconductor OrgDb-sets, the TxDb sets or both. It is obviously important to make sure that the genome build used is the same as used in aligning the RNA-seq data.  

} 
 
\author{Jesper R. Gadin, Lasse Folkersen}
\seealso{
	\itemize{
		\item The \code{\link{ASEset}} class which the locationplot function can be called up on.
	}
}
\examples{

	data(ASEset)
	locationplot(ASEset)

	#SNPs are plotted in the order in which they are found. 
	#This can be sorted according to location as follows:
	locationplot(ASEset[order(start(rowData(ASEset))),])
	
	#for ASEsets with fewer SNPs the 'count' type plot is
	# useful for detailed visualization.
	locationplot(ASEset,type="count",strand="nonStranded")
}
\keyword{locationplot}


