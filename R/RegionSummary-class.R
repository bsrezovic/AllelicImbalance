#'@include initialize-methods.R
NULL

#' RegionSummary class
#' 
#' Object that holds results from AI detection.
#'
#' The RegionSummary-class contains 
#'
#' @name RegionSummary-class
#' @rdname RegionSummary-class
#' @aliases RegionSummary-class RegionSummary RegionSummary-method
#' @docType class
#' @param value argument used for replacement
#' @param ... pass arguments to internal functions
#'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords class RegionSummary
#' @examples
#'
#' #some code
#'
#' @exportClass RegionSummary
NULL

#' @rdname RegionSummary-class
#' @exportClass RegionSummary
setClass("RegionSummary", contains = "RangedSummarizedExperiment",
	representation(
		meta = "list",
		sumnames = "character"
	)
)

#' @rdname ASEset-class
#' @export 
setGeneric("sumnames", function(x, ...){
    standardGeneric("sumnames")
})

#' @rdname ASEset-class
#' @export 
setMethod("sumnames", signature(x = "RegionSummary"), function(x) {
	x@sumnames
})

