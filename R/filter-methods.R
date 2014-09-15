#'@include ASEset-class.R
NULL

#' genotype filter methods 
#' 
#' useful genotype filters
#' 
#' These filters are called upon ASEset objects 
#' 
#' @name genofilters
#' @rdname genofilters
#' @aliases hetFilt hetFilt,ASEset-method 
#' @docType methods
#' @param x \code{ASEset} object
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords filter
#' @examples
#' 
#' #load example data
#' data(ASEset)
#' a <- ASEset
#*
#' genotype(a) <- inferGenotypes(a)
#' hets <- hetFilt(a) 
#' 
#' @exportMethod hetFilt
NULL

# @rdname genofilters
setGeneric("hetFilt", function(x){
    standardGeneric("hetFilt")
})

setMethod("hetFilt", signature(x = "ASEset"), function(x){
	
	 matrix(vapply(
		    genotype(x),
			function(y){substring(y,1,1)==substring(y,3,3)}, numeric(1)),
			nrow=nrow(x), ncol=ncol(x), dimnames=list(rownames(x),colnames(x))
	) == 0

})

