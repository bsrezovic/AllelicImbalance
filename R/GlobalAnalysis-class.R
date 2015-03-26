#'@include initialize-methods.R
NULL

#' GlobalAnalysis class
#' 
#' Object that holds results from a global AI analysis including reference bias
#' estimations and AI detection.
#'
#' The GlobalAnalysis-class contains summaries and "pre-configured and pre-calculated lattice plots" needed to create an AI-report
#'
#' @name GlobalAnalysis-class
#' @rdname GlobalAnalysis-class
#' @aliases GlobalAnalysis-class GlobalAnalysis GlobalAnalysis-method
#' @docType class
#' @param x ASEset object or list of ASEsets
#' @param TxDb A \code{transcriptDb} object
#' @param ... pass arguments to internal functions
#' @return An object of class GlobalAnalysis containing all data to make report.
#'
#'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords class ASEset
#' @examples
#'
#' data(ASEset)
#' #a <- ASEset
#' #gba <- gba(a)
#' 
#' #report(gba)
#' #write.tables(gba)
#' #graphs(gba)
#' #as.list(gba)
#'
#' @exportClass GlobalAnalysis
NULL

#' @rdname GlobalAnalysis-class
#' @exportClass GlobalAnalysis
setClass("GlobalAnalysis",  
	representation(data = "list"))



