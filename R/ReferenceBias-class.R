#'@include initialize-methods.R
NULL

#' ReferenceBias class
#' 
#' Object that holds Reference bias information
#'
#' Reference bias-class contain annotated fractions of the reference allele
#'
#' Shortly the ReferenceBias-class most prominent purpose is quikly to be able to dispatch on
#' methods like, 'plot', 'summary' and similar.
#' 
#' @name ReferenceBias-class
#' @rdname ReferenceBias-class
#' @aliases ReferenceBias-class ReferenceBias ReferenceBias-method
#' @docType class
#' @param x ReferenceBias object
#' @param TxDb A \code{transcriptDb} object
#' @param strand '+','-' or '*'
#' @param return.class class to return
#' @param ... pass arguments to internal functions
#' @return An object of class ReferenceBias storing reference fractions.
#'
#' @section Constructor: refBias(x = ASEset)
#' 
#' \describe{
#' 
#' Arguments: \item{x}{an ASEset object  }
#'
#' }
#'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords class ASEset
#' @examples
#'
#' data(ASEset)
#' a <- ASEset
#' rb <- RBias(a)
#' 
#' @exportClass ReferenceBias
NULL


#' @rdname ReferenceBias-class
#' @exportClass ReferenceBias
setClass("ReferenceBias", contains = "SummarizedExperiment",
	representation(strands = "vector"))

#' @rdname ReferenceBias-class
#' @export
setGeneric("refbiasByAnnotation", function(x,...){
    standardGeneric("refbiasByAnnotation")
})

#' @rdname ReferenceBias-class
#' @export
setMethod("refbiasByAnnotation", signature(x = "ReferenceBias"), function(x,TxDb){
	cat("not implemented yet, sorry")	
})

#hidden function
.selectArrayElement <- function(strand) 
	{
    if (!sum(strand %in% c("+", "-", "*")) > 0) {
        stop("strand parameter has to be either '+', '-', '*' ")
    }
    
    if (strand == "+") {
        el <- 2
    } else if (strand == "-") {
        el <- 3
    } else if (strand == "*") {
        el <- 1
	} else {
        stop("not existing strand option")
    }
	el
}

#' @rdname ReferenceBias-class
#' @export
setMethod("frequency", signature(x = "ReferenceBias"), function(x, 
	return.class = "matrix", strand = "*"
	) {

	el <- .selectArrayElement(strand)
	ar <- assays(x)[["referenceFrequency"]][,,el]

	if(return.class=="matrix"){
		return(ar)
	}else{
		stop("return.class has to be 'matrix' ")
	}
})

#' @rdname ReferenceBias-class
#' @export
setGeneric("hetPerSample", function(x,...){
    standardGeneric("hetPerSample")
})

#' @rdname ReferenceBias-class
#' @export
setMethod("hetPerSample", signature(x = "ReferenceBias"), function(x, 
	return.class = "vector", strand = "*"
	) {

	el <- .selectArrayElement(strand)
	ar <- assays(x)[["referenceFrequency"]][,,el]
	ar[!is.na(ar)] <- 1
	ar[is.na(ar)] <- 0

	#dimnames(ar) <- list(rownames(x),colnames(x),x@strands)

	vec <- apply(ar,2,sum)

	if(return.class=="vector"){
		return(vec)
	}else{
		stop("return.class has to be 'vector' ")
	}

})

#' @rdname ReferenceBias-class
#' @export
setGeneric("hetPerSnp", function(x,...){
    standardGeneric("hetPerSnp")
})

#' @rdname ReferenceBias-class
#' @export
setMethod("hetPerSnp", signature(x = "ReferenceBias"), function(x, 
	return.class = "vector", strand = "*"
	) {

	el <- .selectArrayElement(strand)
	ar <- assays(x)[["referenceFrequency"]][,,el]
	ar[!is.na(ar)] <- 1
	ar[is.na(ar)] <- 0

	#dimnames(ar) <- list(rownames(x),colnames(x),x@strands)

	vec <- apply(ar,1,sum)

	if(return.class=="vector"){
		return(ar)
	}else{
		stop("return.class has to be 'vector' ")
	}
	vec
})


