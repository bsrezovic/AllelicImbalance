#'@include initialize-methods.R
NULL

#' LinkVariantAlmlof class
#' 
#' Object that holds results from AI detection.
#'
#' The LinkVariantAlmlof-class contains 
#'
#' @name LinkVariantAlmlof-class
#' @rdname LinkVariantAlmlof-class
#' @aliases LinkVariantAlmlof-class LinkVariantAlmlof LinkVariantAlmlof-method
#' @docType class
#' @param x LinkVariantAlmlof object
#' @param ... pass arguments to internal functions
#'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords class LinkVariantAlmlof
#' @examples
#'
#' #some code
#'
#' @exportClass LinkVariantAlmlof
NULL

#' @rdname LinkVariantAlmlof-class
#' @exportClass LinkVariantAlmlof
setClass("LinkVariantAlmlof", contains = "RangedSummarizedExperiment",
	representation(
		meta = "list"
	)
)

#' @rdname LinkVariantAlmlof-class
#' @export 
setGeneric("pvalue", function(x, ...){
    standardGeneric("pvalue")
})

#' @rdname LinkVariantAlmlof-class
#' @export 
setMethod("pvalue", signature(x = "LinkVariantAlmlof"), function(x) {
		mcols(x)[["LMCommonParam"]][,8]
})



