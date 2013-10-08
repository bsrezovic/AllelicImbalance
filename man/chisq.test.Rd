\name{chisq.test}
\docType{methods}

\alias{chisq.test,ASEset-method}

\title{chi-square test}
\description{
	Performs a chisq.test on an ASEset object.
}
\usage{
#snpAFTotDf comes from the mergeAFList function
\S4method{chisq.test}{ASEset}(x,y="nonStranded")
}
\arguments{
	\item{x}{ \code{ASEset} object}
	\item{y}{strand option}
}
\details{
	The test is performed on one strand in an ASEset object.
}
\value{
        \code{chisq.test} returns a matrix with the chisq.test P-value for each SNP and sample
}
\author{Jesper R. Gadin, Lasse Folkersen}

\seealso{
	\itemize{
		\item The \code{\link{binom.test}} which is another test that can be applied on an \link{ASEset} object.
	}
}

\examples{
	#load example data
	data(ASEset)
	
	#make a chi-square test on default nonStranded strand 
	chisq.test(ASEset)
} 	

\keyword{ chi-square test}

