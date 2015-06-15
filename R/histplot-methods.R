#' histogram plots
#' 
#' uses base graphics hist plot
#' 
#' The histogram will show the density over frequencies for each sample
#' 
#' @name histplot
#' @rdname histplot
#' @aliases hist hist,ReferenceBias-method 
#' @docType methods
#' @param x \code{ReferenceBias} object or \code{ASEset} object
#' @param strand '+','-' or '*'
#' @param type 'mean' (only one option atm)
#' @param log an integer to log each value (integer 10 for log10)
#' @param ... arguments to forward to interal boxplots function
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords plot hist
#' @examples
#' 
#' ##load example data
#' 
#' #data(ASEset)
#' #a <- ASEset
#' #hist(a)
#' 
NULL

#' @rdname histplot
setMethod("hist", signature(x = "ASEset"), function(x, strand="*", type="mean", log=1, ...){

		if(log==1){
			counts <- countsPerSnp(x, strand=strand, return.type=type, return.class="vector")
		}else{
			counts <- log(countsPerSnp(x, strand=strand, return.type=type, return.class="vector"),log)
		}

		hi <- hist(counts,breaks = 40, freq=TRUE, ...)
		invisible(hi)
})


