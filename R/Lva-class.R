#'@include initialize-methods.R
NULL

#' Lva class
#' 
#' Object that holds results from AI detection.
#'
#' The Lva-class contains 
#'
#' @name Lva-class
#' @rdname Lva-class
#' @aliases Lva-class Lva Lva-method
#' @docType class
#' @param value argument used for replacement
#' @param ... pass arguments to internal functions
#'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords class Lva
#' @examples
#'
#' #some code
#'
#' @exportClass Lva
NULL

#' @rdname Lva-class
#' @exportClass Lva
setClass("Lva", contains = "RangedSummarizedExperiment",
	representation(
		meta = "list"
	)
)



