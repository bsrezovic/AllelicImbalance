\name{ASEset-lbarplot}
\docType{methods}
\alias{ASEset-lbarplot}
\alias{lbarplot}
\alias{lbarplot,ASEset-method}


\title{lbarplot ASEset objects}
\description{
	Generates lbarplots for ASEset objects. Two levels of plotting detail are provided: a detailed lbarplot of read counts by allele useful for fewer samples and SNPs, and a less detailed lbarplot of the fraction of imbalance, useful for more samples and SNPs. 
}
\usage{
  \S4method{lbarplot}{ASEset}(
	x,
	type="count",
	strand="+",
	mainVec=rep("",nrow(x)),
	verbose=FALSE,
	...)
}
\arguments{
	\item{x}{An \code{ASEset} object}
	\item{type}{"count" or "fraction"}
	\item{strand}{Five options, "nonStranded","+", "-", "both" or "*"}
	\item{mainVec}{text to use as main label}
	\item{verbose}{Makes function more talkative}
	\item{...}{for simpler generics when extending function}
}	


\details{
	This function serves the same purpose as the normal barplot, but with trellis graphics using lattice, to be able to integrate well with Gviz track functionality.	
}

\author{Jesper R. Gadin}

\seealso{
	\itemize{
		\item The \code{\link{ASEset}} class which the lbarplot function can be called up on.
		\item The \code{\link{barplot}} non trellis barplot.
	}
}

\examples{
	data(ASEset)
	lbarplot(ASEset[1])
}
\keyword{lbarplot}

