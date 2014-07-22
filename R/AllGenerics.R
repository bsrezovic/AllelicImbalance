#accessors
setGeneric("alleleCounts", function(x,strand="nonStranded") {standardGeneric("alleleCounts")})
setGeneric("mapBias", function(x) {standardGeneric("mapBias")})
setGeneric("fraction", function(x, strand="nonStranded", verbose=FALSE) {standardGeneric("fraction")})
setGeneric("arank", function(x, ret="names", strand="nonStranded", ... ) {standardGeneric("arank")})

#visuals
setGeneric("barplot")
setGeneric("glocationplot", function(x,
		type="fraction",
		strand="nonStranded")
		{standardGeneric("glocationplot")})
				     	
setGeneric("locationplot", function(x,
		type="fraction",
		strand="nonStranded",
		yaxis=TRUE,
		xaxis=FALSE,
		xlab=FALSE,
		ylab=TRUE,
		legend.colnames = "", 
		size=0.9,
		main=NULL,
		pValue=FALSE,
		cex.main=0.7,
		cex.ylab=0.6,
		cex.legend= 0.6,
		OrgDb=NULL,
		TxDb=NULL,
		verbose=TRUE,
		...) {standardGeneric("locationplot")})

#tests
setGeneric("chisq.test")
setGeneric("binom.test")





