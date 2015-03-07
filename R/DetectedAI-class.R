#'@include initialize-methods.R
NULL

#' DetectedAI class
#' 
#' Object that holds results from AI detection.
#'
#' The DetectedAI-class contains 
#'
#' @name DetectedAI-class
#' @rdname DetectedAI-class
#' @aliases DetectedAI-class DetectedAI DetectedAI-method
#' @docType class
#' @param x ASEset object or list of ASEsets
#' @param return.class type of class returned eg. "list or ""array".
#' @param ... pass arguments to internal functions
#' @return An object of class DetectedAI containing logics for thresholds of interest.
#'
#'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords class ASEset
#' @examples
#'
#' data(ASEset)
#' a <- ASEset
#' dai <- detectAI(a)
#' 
#' #summary(gba)
#' #write.tables(dai)
#'
#' @exportClass DetectedAI
NULL

#' @rdname DetectedAI-class
#' @exportClass DetectedAI
setClass("DetectedAI", contains = "SummarizedExperiment",
	representation(
		strand = "character",
		threshold.count.sample.names="character",
		threshold.frequency.names="character",
		threshold.delta.frequency.names="character",
		threshold.pvalue.names="character"
	)
)

#' @rdname DetectedAI-class
#' @export 
setGeneric("referenceFrequency", function(x, ...) {
    standardGeneric("referenceFrequency")
})

#' @rdname DetectedAI-class
#' @export 
setMethod("referenceFrequency", signature(x = "DetectedAI"), function(x,
	return.class="array") {

	if(return.class=="array"){
		return(assays(x)[["reference.frequency"]])
	}
})

#' @rdname DetectedAI-class
#' @export 
setGeneric("thresholdFrequency", function(x, ...) {
    standardGeneric("thresholdFrequency")
})

#' @rdname DetectedAI-class
#' @export 
setMethod("thresholdFrequency", signature(x = "DetectedAI"), function(x,
	return.class="array") {

	if(return.class=="array"){
		return(assays(x)[["threshold.frequency"]])
	}
})

#' @rdname DetectedAI-class
#' @export 
setGeneric("thresholdCountSample", function(x, ...) {
    standardGeneric("thresholdCountSample")
})

#' @rdname DetectedAI-class
#' @export 
setMethod("thresholdCountSample", signature(x = "DetectedAI"), function(x,
	return.class="array") {

	if(return.class=="array"){
		return(assays(x)[["threshold.count.sample"]])
	}
})

#' @rdname DetectedAI-class
#' @export 
setGeneric("thresholdDeltaFrequency", function(x, ...) {
    standardGeneric("thresholdDeltaFrequency")
})

#' @rdname DetectedAI-class
#' @export 
setMethod("thresholdDeltaFrequency", signature(x = "DetectedAI"), function(x,
	return.class="array") {

	if(return.class=="array"){
		return(assays(x)[["threshold.delta.frequency"]])
	}
})

#' @rdname DetectedAI-class
#' @export 
setGeneric("thresholdPvalue", function(x, ...) {
    standardGeneric("thresholdPvalue")
})

#' @rdname DetectedAI-class
#' @export 
setMethod("thresholdPvalue", signature(x = "DetectedAI"), function(x,
	return.class="array") {

	if(return.class=="array"){
		return(assays(x)[["threshold.pvalue"]])
	}
})

