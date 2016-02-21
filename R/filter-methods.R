#'@include ASEset-class.R
NULL

#' genotype filter methods 
#' 
#' useful genotype filters
#' 
#' hetFilt returns TRUE if the samples is heterozygote, based on stored genotype information
#' present in the phase data.
#' 
#' @name ASEset-filters
#' @rdname ASEset-filters
#' @aliases hetFilt hetFilt,ASEset-method 
#' @docType methods
#' @param x ASEset object
#' @param source 'genotype' or 'alleleCounts'
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
NULL

#' @rdname ASEset-filters
#' @export 
setGeneric("hetFilt", function(x, ...){
    standardGeneric("hetFilt")
})

#' @rdname ASEset-filters
#' @export 
setMethod("hetFilt", signature(x = "ASEset"), 
	function(x, source="genotype", ...)
	{
		if(source=="genotype"){
			.heterozygozityFromPhaseArray(phase(x , return.class="array"))
		}else if(source=="alleleCounts"){
			stop("not implemented")
		}else{
			stop("source must be 'genotype' or 'alleleCounts' ")
		}
	}
)

### -------------------------------------------------------------------------
### helpers for hetFilt
###
.heterozygozityFromPhaseArray <- function(x){
	!(x[,,1] == x[,,2])
}

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
NULL

#' @rdname ASEset-filters
#' @export 
setGeneric("multiAllelicFilt", function(x, ...){
    standardGeneric("multiAllelicFilt")
})

#' @rdname ASEset-filters
#' @export 
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

#' @rdname ASEset-filters
#' @export 
setGeneric("minFreqFilt", function(x, ...){
    standardGeneric("minFreqFilt")
})

#' @rdname ASEset-filters
#' @export 
setMethod("minFreqFilt", signature(x = "ASEset"), 
		function(x, strand="*", threshold.frequency=0.10, replace.with="zero", 
				 return.class="ASEset", sum="all" )
				
		{
		thr <- threshold.frequency
		fr <- frequency(x, strand=strand, threshold.count.sample=1,
				return.class="array")

		#create tfmat
		if(sum=="each"){
			#check if ref and alt allele are present
			if(!refExist(x)){stop("ref(x) is empty")}
			if(!altExist(x)){stop("alt(x) is empty")}
			#check if ref and alt are character, otherwise give a warning and 
			#coerce to character
			ref <- .verboseCoerceToCharacter(ref(x)) 
			alt <- .verboseCoerceToCharacter(alt(x)) 
			tfmat <- .toKeepMatrixMinFreqFilterEach(fr,ref,alt,ncol(x),thr,x@variants)

		}else if(sum=="all"){
			tfmat <- .Na2False(.toKeepMatrixMinFreqFilterAll(fr, thr))
		}

		#return object
		if(return.class=="ASEset"){
			#will only filter out the selected strand, and while the unknown strand is
			#based on the plus and minus we have to filter all of them, to have the back
			#free. In the future the assays(x)[["countsUnknown"]] will be dropped so it 
			#always calculated on the fly as the sum of the plus and minus strand
			tfarray <- .expandMatrixToArray(tfmat, length(x@variants))
			if(strand=="*"){
				if(.unknownStrandCountsExists(x)) alleleCounts(x, strand="*", return.class="array")[!tfarray] <- 0
				if(.minusStrandCountsExists(x)) alleleCounts(x, strand="-", return.class="array")[!tfarray] <- 0
				if(.plusStrandCountsExists(x)) alleleCounts(x, strand="+", return.class="array")[!tfarray] <- 0
			}else if(strand=="+" | strand=="-"){
				alleleCounts(x, strand=strand, return.class="array")[!tfarray] <- 0
			}
			return(x)
		}
		if(return.class=="array") return(.expandMatrixToArray(tfmat, length(x@variants)))
		if(return.class=="matrix") return(tfmat)
})
### -------------------------------------------------------------------------
### helpers for minFreqFilt
###

#get the replacement matrix for "all"
.toKeepMatrixMinFreqFilterAll <- function(fr, thr){
		tf <- fr >= thr
		apply(tf, c(1,2), function(x){sum(x)>=2})	
}

#get the replacement matrix for "each"
.toKeepMatrixMinFreqFilterEach <- function(fr,ref,alt,nc,thr,var){
			#check that overall dimensions are fine
			if(!length(nc)==1) stop("dimensions do not agree in .altAlleleFreqThreshold")
			if(any(!dim(fr) == c(length(ref), nc, length(var))) | 
			   !(length(ref) == length(alt)) |
			   !(length(nc) == 1) |
			   !(length(thr) == 1)
			) stop("dimensions do not agree in .altAlleleFreqThreshold")
			

			#take out count values only for ref allele
			ref.ar <- .arrayFromAlleleVector(var, ref, nc )
			alt.ar <- .arrayFromAlleleVector(var, alt, nc )

			#use on counts
			ref.fr <- .subsetArrayToMatrix(fr, ref.ar)
			alt.fr <- .subsetArrayToMatrix(fr, alt.ar)

			#at least eg. 0.1 in each group
			.Na2False(ref.fr > thr) & .Na2False(alt.fr > thr)

			#tfFilt <- .Na2False(ref.fr > thr) & .Na2False(alt.fr > thr)
			#take out these pairs and look at them (manually control that all is fine)
			#this is atm not tested in a unit test, but should be covered by
			#another test
			#matrix(c(ref.fr[tfFilt], alt.fr[tfFilt]), ncol=2)
			#matrix(c(ref.fr[!tfFilt], alt.fr[!tfFilt]), ncol=2)

			#extend the matrix to 3d for the subset/replacement now done outside this unit
			#.expandMatrixToArray(tfFilt, length(var))
}


#' @rdname ASEset-filters
#' @export 
setGeneric("minCountFilt", function(x, ...){
    standardGeneric("minCountFilt")
})

#' @rdname ASEset-filters
#' @export 
setMethod("minCountFilt", signature(x = "ASEset"), 
		function(x, strand="*", threshold.counts=1,
				   sum="all", replace.with="zero", return.class="ASEset"){

		#set shorter name
		thr <- threshold.counts

		#extract alleleCounts
		ac <- alleleCounts(x, strand=strand, return.class="array")

		#create tfmat
		if(sum=="all") tfmat <- .toKeepMatrixMinCountFilterAll(ac, thr)
		else if(sum=="each"){
			#check if ref and alt allele are present
			if(!refExist(x)){stop("ref(x) is empty")}
			if(!altExist(x)){stop("alt(x) is empty")}
			#check if ref and alt are character, otherwise give a warning and 
			#coerce to character
			ref <- .verboseCoerceToCharacter(ref(x)) 
			alt <- .verboseCoerceToCharacter(alt(x)) 
			tfmat <- .toKeepMatrixMinCountFilterEach(ac,ref,alt,ncol(x),thr,x@variants)
		}

		#return object
		if(return.class=="ASEset"){
			#will only filter out the selected strand, and while the unknown strand is
			#based on the plus and minus we have to filter all of them, to have the back
			#free. In the future the assays(x)[["countsUnknown"]] will be dropped so it 
			#always calculated on the fly as the sum of the plus and minus strand
			tfarray <- .expandMatrixToArray(tfmat, length(x@variants))
			if(strand=="*"){
				if(.unknownStrandCountsExists(x)) alleleCounts(x, strand="*", return.class="array")[!tfarray] <- 0
				if(.minusStrandCountsExists(x)) alleleCounts(x, strand="-", return.class="array")[!tfarray] <- 0
				if(.plusStrandCountsExists(x)) alleleCounts(x, strand="+", return.class="array")[!tfarray] <- 0
			}else if(strand=="+" | strand=="-"){
				alleleCounts(x, strand=strand, return.class="array")[!tfarray] <- 0
			}
			return(x)
		}
		if(return.class=="array") return(.expandMatrixToArray(tfmat, length(x@variants)))
		if(return.class=="matrix") return(tfmat)

})

### -------------------------------------------------------------------------
### helpers for minCountFilt
###

#get the replacement matrix for "all"
.toKeepMatrixMinCountFilterAll <- function(ac, thr){
	apply(ac, c(1,2), function(x,thr){sum(x)> thr},thr=thr)
}

#get the replacement matrix for "each"
.toKeepMatrixMinCountFilterEach <- function(ac,ref,alt,nc,thr,var){
			#check that overall dimensions are fine
			if(!length(nc)==1) stop("dimensions do not agree in .toKeepArrayMinCountFilter")
			if(any(!dim(ac) == c(length(ref), nc, length(var))) | 
			   !(length(ref) == length(alt)) |
			   !(length(nc) == 1) |
			   !(length(thr) == 1)
			) stop("dimensions do not agree in .toKeepArrayMinCountFilter")
			

			#take out count values only for ref allele
			ref.ar <- .arrayFromAlleleVector(var, ref, nc )
			alt.ar <- .arrayFromAlleleVector(var, alt, nc )

			#use on counts
			ref.ac <- .subsetArrayToMatrix(ac, ref.ar)
			alt.ac <- .subsetArrayToMatrix(ac, alt.ar)

			#at least 30 in both groups
			(ref.ac > thr) & (alt.ac > thr)

			#take out these pairs and look at them (manually control that all is fine)
			#this is atm not tested in a unit test, but should be covered by
			#another test
			#matrix(c(ref.ac[tfFilt], alt.ac[tfFilt]), ncol=2)
			#matrix(c(ref.ac[!tfFilt], alt.ac[!tfFilt]), ncol=2)

			#extend the matrix to 3d for the subset/replacement now done outside this unit
			#.expandMatrixToArray(tfFilt, length(var))
}

