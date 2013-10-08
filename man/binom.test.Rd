\name{binom.test}
\docType{methods}

\alias{binom.test,ASEset-method}

\title{binomial test}
\description{
	Performs a binomial test on an ASEset object.
}
\usage{
#snpAFTotDf comes from the mergeAFList function
\S4method{binom.test}{ASEset}(x,n="nonStranded")
}
\arguments{
	\item{x}{ \code{ASEset} object}
	\item{n}{strand option}
}
\details{
	the test can only be applied to one strand at the time.
}
\value{
        \code{binom.test}  returns a matrix 
}
\author{Jesper R. Gadin, Lasse Folkersen}


\seealso{
	\itemize{
		\item The \code{\link{chisq.test}} which is another test that can be applied on an \link{ASEset} object.
	}
}

\examples{
	#load example data
	data(ASEset)
	
	#make a binomial test
	binom.test(ASEset,"nonStranded")

}

\keyword{ binomial test }

