#' DetectedAI summary
#' 
#' Summary helper functions for the DetectedAI-class
#' 
#' Summary helper functions. The documentation will
#' be improved before next release.
#' 
#' @name DetectedAI-summary
#' @rdname DetectedAI-summary
#' @aliases frequency_vs_threshold_variable_summary frequency_vs_threshold_variable_summary,DetectedAI-method detectedAI_vs_threshold_variable_summary detectedAI_vs_threshold_variable_summary,DetectedAI-method
#' @param x detectedAI object 
#' @param var string, see details for available options
#' @param ... pass on variables internally
#' @param return.class 'matrix' or 'array'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords list
#' @examples
#' 
#' #some example code here
#' #generate example
#' data(ASEset)
#' a <- ASEset
#' dai <- detectAI(a, 
#' 			threshold.count.sample=1:50,
#' 			threshold.frequency=seq(0,0.5,by=0.01),
#' 			threshold.delta.frequency=seq(0,0.5,by=0.01),
#' 			threshold.pvalue=rev(seq(0.001,0.05, by=0.005))
#' )
#' 
#' frequency_vs_threshold_variable_summary(dai)
#' 
#' 
NULL

#' @rdname DetectedAI-summary
#' @export
setGeneric("frequency_vs_threshold_variable_summary", function(x, ... 
	){
    standardGeneric("frequency_vs_threshold_variable_summary")
})

#' @rdname DetectedAI-summary
#' @export
setMethod("frequency_vs_threshold_variable_summary", signature(x = "DetectedAI"),
		function(x,var="threshold.count.sample", return.class="matrix", ...){

		#check if assay is present
		fr <- assays(x)[["reference.frequency"]]
		ar.var <- assays(x)[[var]]
		ar.fr <- array(fr,dim=c(nrow(fr),ncol(fr),dim(ar.var)[3]),
					  dimnames=list(rownames(x),colnames(x),NULL) )

		is.na(ar.var) <- FALSE 
		ar.fr[!ar.var] <- NA

		#colSums(ar.fr,na.rm=TRUE)

		if(return.class=="matrix"){
			apply(ar.fr,c(2, 3),mean,na.rm=TRUE)
		}else if(return.class=="array"){
			ar.fr	
		}

})

#' @rdname DetectedAI-summary
#' @export
setGeneric("detectedAI_vs_threshold_variable_summary", function(x, ... 
	){
    standardGeneric("detectedAI_vs_threshold_variable_summary")
})

#' @rdname DetectedAI-summary
#' @export
setMethod("detectedAI_vs_threshold_variable_summary", signature(x = "DetectedAI"),
		function(x, var="threshold.count.sample"){

	#check if assay is present
	apply(assays(x)[[var]],c(2, 3),sum,na.rm=TRUE)

})

#' @rdname DetectedAI-summary
#' @export
setGeneric("usedSNPs_vs_threshold_variable_summary", function(x, ... 
	){
    standardGeneric("usedSNPs_vs_threshold_variable_summary")
})

#' @rdname DetectedAI-summary
#' @export
setMethod("usedSNPs_vs_threshold_variable_summary", signature(x = "DetectedAI"),
		function(x, var="threshold.count.sample"){

	# filter variable for heterozygotes used in calculation
	# based on non NA:s in the reference.frequency assay and threshold variable

	tf <- array(!is.na(assays(x)[["reference.frequency"]]) & !(assays(x)[["reference.frequency"]]==1) & !(assays(x)[["reference.frequency"]]==0), dim=dim(assays(x)[[var]])) & assays(x)[[var]]

	apply(tf,c(2, 3),sum,na.rm=TRUE)

})


