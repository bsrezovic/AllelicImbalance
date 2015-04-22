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
#'
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

#' multi-allelic filter methods 
#' 
#' filter on multiallelic snps
#' 
#' based on the allele counts for all four variants A, T, G and C and returns true 
#' if there is counts enough suggesting a third or more alleles. The sensitivity can 
#' be specified using 'threshold.count.sample' and 'threshold.frequency'.
#' 
#' @name multiAllelicFilt
#' @rdname genofilters
#' @aliases multiAllelicFilt multiAllelicFilt,ASEset-method 
#' @docType methods
#' @param x \code{ASEset} object
#' @param strand strand to infer from
#' @param inferOver 'eachSample' or 'allSamples' 
#' @param threshold.count.sample least amount of counts to try to infer allele
#' @param threshold.frequency least fraction to classify (see details)
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords filter
#' @examples
#' 
#' #load example data
#' data(ASEset)
#' a <- ASEset
#'
#' multiAllelicFilt(a)
#' 
#' @exportMethod multiAllelicFilt
NULL

# @rdname multiAllelicFilt
setGeneric("multiAllelicFilt", function(x, ...){
    standardGeneric("multiAllelicFilt")
})

setMethod("multiAllelicFilt", signature(x = "ASEset"), 
		function(x, strand="*", threshold.count.sample=10, threshold.frequency=0.10,
				   filterOver="eachSample"){

		#find SNP types
		fr <- frequency(x, strand=strand, 
			threshold.count.sample=threshold.count.sample,
			return.class="array")

		#apply over alleles
		alleles <- apply(fr,c(1,2), function(x){
			sum(x >= threshold.frequency)
		})

		if(filterOver=="allSamples"){
			apply(alleles,1, function(x){
				any(x > 2 ,na.rm=TRUE)
			})
		}else if(filterOver=="eachSample"){
			alleles > 2
		}
})


