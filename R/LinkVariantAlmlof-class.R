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
#' @param value argument used for replacement
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



