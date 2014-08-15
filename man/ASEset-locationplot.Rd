\name{ASEset-glocationplot}
\docType{methods}
\alias{ASEset-glocationplot}
\alias{glocationplot}
\alias{glocationplot,ASEset-method}


\title{glocationplot ASEset objects}
\description{
	plotting ASE effects over a specific genomic region using Gviz functionality
}
\usage{
\S4method{glocationplot}{ASEset}(
	x,
		type="fraction",
		strand="nonStranded",
		BamGAL=NULL,
		...
	)
}
\arguments{
	\item{x}{ an ASEset object.}
	\item{type}{"fraction" or "count"}
	\item{strand}{"+","-","both". This argument determines which strand is plotted. See \code{getAlleleCounts} for more information on strand.}
	\item{...}{arguments passed on to barplot function}
}	

\details{
	The glocationplot methods visualises the distribution of ASE over a larger region on one chromosome. It takes and ASEset object as well as additional information on plot type (see \code{\link{lbarplot}}), strand type (see \code{\link{getAlleleCounts}}), Annotation tracks are created from the Gviz packageh. It is obviously important to make sure that the genome build used is set correctly, e.g. "hg19".

} 
 
\author{Jesper R. Gadin}
\seealso{
	\itemize{
		\item The \code{\link{ASEset}} class which the glocationplot function can be called up on.
	}
}
\examples{
	data(ASEset)
	glocationplot(ASEset,strand="+")
	
	#for ASEsets with fewer SNPs the 'count' type plot is useful 
	glocationplot(ASEset,type="count",strand="+")
}
\keyword{glocationplot}


