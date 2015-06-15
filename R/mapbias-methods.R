#' Reference allele
#' 
#' Extract the allele based on SNP location from the reference fasta file
#' 
#' The alleles will be placed in the rowRanges() meta column 'ref'
#' 
#' 
#' @name refAllele 
#' @rdname refAllele
#' @aliases refAllele,ASEset-method
#' @docType methods
#' @param x \code{ASEset} object
#' @param fasta path to fasta file, index should be located in the same folder
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords reference mapbias
#' @examples
#' 
#' #load example data
#' data(ASEset.sim)
#'
#' fasta <- system.file('extdata/hg19.chr17.subset.fa', package='AllelicImbalance')
#' a <- refAllele(ASEset.sim,fasta=fasta) 
#'
#' @exportMethod refAllele
NULL

#' @rdname refAllele
setGeneric("refAllele", function(x, fasta ){
    standardGeneric("refAllele")
})

setMethod("refAllele", signature(x = "ASEset"), function(x, fasta){

	#does the fasta file exist?
	if(!file.exists(fasta)){
		stop("fasta file doesnt exist")
	}

	#make FaFile
	fl <- FaFile(fasta)

	#check if the index file is present otherwise tell the user to use the
	#indexFa(FaFile("pathToReference")) command
	if(!file.exists(paste(fasta,".fai",sep=""))){
		cat("could not find index file\n")
		cat("creates a new index file")
		indexFa(FaFile(fasta))
		cat("finished creating new index file")
	}
	#IMPORTANT! The  index command only needs to be executed once
	#indexFa(fl) #creates a new file as index in the same directory but
	#with extension *.fai 

	#open,scan,close file
	open(fl)
	#check if seqleveles are present
	fa.info <- scanFaIndex(fl)
	if(!all(seqlevels(x) %in% seqlevels(fa.info))){
		stop("seqlevels in object x are not in fasta index file")
	}
	
	ref <- scanFa(fl, param=rowRanges(x))
	close(fl)
	
	#mcols(x)[["ref"]] <- as.vector(ref)
	#x
	as.vector(ref)
})

#' Generate default mapbias from genotype
#' 
#' Create mapbias array from genotype matrix requires genotype information
#' 
#' Default mapbias will be 0.5 for bi-allelic snps and 1 for
#' homozygots. For genotypes with NA, 0.5 will be placed on all four alleles.
#' Therefore tri-allelic can not be used atm. Genotype information has to be 
#' placed in the genotype(x) assay.
#' 
#' @name defaultMapBias
#' @rdname defaultMapBias
#' @aliases defaultMapBias,ASEset-method
#' @docType methods
#' @param x \code{ASEset} object
#' @param return.class "array" or "ASEset"
#' @param ... internal arguments
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords mapbias
#' @examples
#' 
#' #load example data
#' data(ASEset.sim)
#'
#' fasta <- system.file('extdata/hg19.chr17.subset.fa', package='AllelicImbalance')
#' refAllele(ASEset.sim,fasta=fasta)
#' a <- refAllele(ASEset.sim,fasta=fasta) 
#'
NULL

#' @rdname defaultMapBias
#' @export
setGeneric("defaultMapBias", function(x,... ){
    standardGeneric("defaultMapBias")
})

#' @rdname defaultMapBias
#' @export
setMethod("defaultMapBias", signature(x = "ASEset"), function(x, return.class="array"){

		if (!("genotype" %in% names(assays(x)))) {
			stop(paste("genotype matrix is not present as assay in",
					   " ASEset object, see '?inferGenotypes' "))
		}

		g <- genotype(x)
		allele1 <- sub("/.*","",as.character(g))
		allele2 <- sub(".*/","",as.character(g))

		#set NA allele as "A"
		allele1.nona <-  allele1
		allele1.nona[is.na(allele1.nona)] <- "A"
		allele2.nona <-  allele2
		allele2.nona[is.na(allele2.nona)] <- "A"

		ar <- array(0, dim=c(nrow(x),ncol(x),4),
			dimnames=list(rownames(x),colnames(x),x@variants ))

		ar[,,"A"][allele1.nona=="A"] <- ar[,,"A"][allele1.nona=="A"] + 0.5
		ar[,,"C"][allele1.nona=="C"] <- ar[,,"C"][allele1.nona=="C"] + 0.5
		ar[,,"T"][allele1.nona=="T"] <- ar[,,"T"][allele1.nona=="T"] + 0.5
		ar[,,"G"][allele1.nona=="G"] <- ar[,,"G"][allele1.nona=="G"] + 0.5

		ar[,,"A"][allele2.nona=="A"] <- ar[,,"A"][allele2.nona=="A"] + 0.5
		ar[,,"C"][allele2.nona=="C"] <- ar[,,"C"][allele2.nona=="C"] + 0.5
		ar[,,"T"][allele2.nona=="T"] <- ar[,,"T"][allele2.nona=="T"] + 0.5
		ar[,,"G"][allele2.nona=="G"] <- ar[,,"G"][allele2.nona=="G"] + 0.5

		if(return.class=="array"){
			ar
		}else if(return.class=="ASEset"){
			assays(x)[["mapBias"]] <- ar
			x
		}
})

#' Random ref allele from genotype
#' 
#' Create a vector of random reference alleles
#' 
#' Randomly shuffles which of the two alleles for each genotype that is 
#' indicated as reference allele, based on either allele count information
#' or previous ref and alt alleles.
#' 
#' When the source is 'alleleCounts', the two most expressed alleles are taken
#' as reference and alternative allele. 
#' 
#' @name randomRef
#' @rdname randomRef
#' @aliases randomRef,ASEset-method
#' @docType methods
#' @param x \code{ASEset} object
#' @param ... internal arguments
#' @param source 'alleleCounts'
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords mapbias
#' @examples
#' 
#' #load example data
#' data(ASEset.sim)
#' a <- ASEset.sim
#'
#' ref(a) <- randomRef(a, source = 'alleleCounts') 
#'
NULL

#' @rdname randomRef
#' @export
setGeneric("randomRef", function(x,... ){
    standardGeneric("randomRef")
})

#' @rdname randomRef
#' @export
setMethod("randomRef", signature(x = "ASEset"), 
	function(x, source="alleleCounts", ...)
{
	#check for old arguments
	if("inferGenotypes" %in% list(...)){
		stop("inferGenotypes is not an available argument anymore")
	}

	if(source=="alleleCounts"){
		apply(arank(x, return.class="matrix")[,1:2],1,sample,1)
#	}else if(source=="refAndAlt"){
#		
#		if(!("ref" %in% names(mcols(x)))){
#			stop(paste("ref allele is not present in mcols in",
#					   " ASEset object, see '?ASEset' "))
#		}
#		if(!("alt" %in% names(mcols(x)))){
#			stop(paste("alt allele is not present in mcols in",
#					   " ASEset object, see '?inferAltAllele' "))
#		}
#
#		matrix(c(ref(x),alt(x)), nrow(x),ncol=2)
#
	}else{
		stop("selected source is not available")
	}
})

#' mapBias for reference allele 
#' 
#' Create a matrix of bias for the reference allele
#' 
#' select the expected frequency for the reference allele
#' 
#' @name mapBiasRef
#' @rdname mapBiasRef
#' @aliases mapBiasRef,ASEset-method
#' @docType methods
#' @param x \code{ASEset} object
#' @param ... internal arguments
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords mapbias
#' @examples
#' 
#' #load example data
#' data(ASEset)
#' a <- ASEset
#'
#' mat <- mapBiasRef(a) 
#'
NULL

#' @rdname mapBiasRef
#' @export
setGeneric("mapBiasRef", function(x,... ){
    standardGeneric("mapBiasRef")
})

#' @rdname mapBiasRef
#' @export
setMethod("mapBiasRef", signature(x = "ASEset"), function(x){

	#check presence of mapBias array
    if (!("mapBias" %in% names(assays(x)))) {
		stop("column name 'ref' in mcols(x) is required")
		stop(paste("genotype matrix is not present as assay in",
				   " ASEset object, see '?inferGenotypes' ",
				  )) 
	}

	#check presence of reference allele
	if(!("ref" %in% colnames(mcols(x)))){
		stop("column name 'ref' in mcols(x) is required")
	}

	matrix(aperm(mapBias(x, return.class="array"),c(3,2,1))[aperm(array(matrix(
		x@variants, ncol=length(x@variants),
		 nrow=nrow(x), byrow=TRUE) == mcols(x)[,"ref"]
	 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
	],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))

})


