#'@include initialize-methods.R
NULL

#' RiskVariant class
#' 
#' Object that holds results from AI detection.
#'
#' The RiskVariant-class contains 
#'
#' @name RiskVariant-class
#' @rdname RiskVariant-class
#' @aliases RiskVariant-class RiskVariant RiskVariant-method
#' @docType class
#' @param x RiskVariant object or list of RiskVariants
#' @param return.class type of class returned eg. "list or ""array".
#' @param value argument used for replacement
#' @param ... pass arguments to internal functions
#'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords class RiskVariant
#' @examples
#'
#' #some code
#'
#' @exportClass RiskVariant
NULL

#' @rdname RiskVariant-class
#' @exportClass RiskVariant
setClass("RiskVariant", contains = "RangedSummarizedExperiment",
	representation(
		meta = "list"
	)
)

#' @rdname RiskVariant-class
#' @export 
setMethod("ref", signature(x = "RiskVariant"), function(x) {

		mcols(x)[["ref"]]
	
})

#' @rdname RiskVariant-class
#' @export 
setMethod("ref<-", signature(x = "RiskVariant"), function(x, value) {

	#we want to accept more classes than just character (e.g. characterList or DNAStringSetList)
	#if(class(value)=="character") {
		#mcols(x)[["ref"]] <- value
	#}else{
		#stop("wrong class")
	#}
	
	mcols(x)[["ref"]] <- value
	x
})


#' @rdname RiskVariant-class
#' @export 
setMethod("alt", signature(x = "RiskVariant"), function(x) {

		mcols(x)[["alt"]]
	
})

#' @rdname RiskVariant-class
#' @export 
#could be renamed to countsAllAlleles
setMethod("alt<-", signature(x = "RiskVariant"), function(x, value) {

	#we want to accept more classes than just character (e.g. characterList or DNAStringSetList)
	#if(class(value)=="character") {

	#	mcols(x)[["alt"]] <- value
	#}else{

	#	stop("wrong class")
	#}
	mcols(x)[["alt"]] <- value
	
	x
})

#' @rdname RiskVariant-class
#' @export 
setMethod("phase", signature(x = "RiskVariant"), function(x, 
	return.class = "matrix" ) {

	if(return.class=="matrix"){
		mat <- phaseArray2phaseMatrix(assays(x)[["phase"]])
		colnames(mat) <- colnames(x)
		rownames(mat) <- rownames(x)
		mat
	}else if(return.class=="array"){
		assays(x)[["phase"]] 
	}
})

#' @rdname RiskVariant-class
#' @export 
setMethod("phase<-", signature(x = "RiskVariant"), function(x,value) {

	if(class(value)=="matrix") {

		if(!identical(dim(x),dim(value))){
			stop("dimension of value does not correspond to the values of object RiskVariant")	
		}
		
		assays(x)[["phase"]] <- phaseMatrix2Array(value, dimnames=NULL)

	}else if(class(value)=="array"){
		assays(x)[["phase"]] <- value
	}
	
	x
})

#####################
# addColnames might not be needed when phase is required

##' @rdname ASEset-class
##' @export 
#setGeneric("addColnames<-", function(x, value){
#    standardGeneric("addColnames<-")
#})
#
##' @rdname RiskVariant-class
##' @export 
#setMethod("addColnames<-", signature(x = "RiskVariant"), function(x,value) {
#
#		  if(is.null(colnames(x))){
#				sset <- SummarizedExperiment(
#					assays = assays(x), 
#					rowRanges = granges(x),
#					colData=DataFrame(row.names=value))
#				rownames(sset) <- names(x)
#				cnrv <- new("RiskVariant", sset,
#					meta = list()
#				)
#
#		  }else{
#			  (stop("object has already colnames, use cbind."))
#		  }
#})
#
