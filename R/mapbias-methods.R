#' reference fraction
#' 
#' The fractions for all heterozygote reference alleles
#' 
#' Neccessary to measure the effect of mapbias over heterozygous SNPs
#'
#' @name refFraction
#' @rdname refFraction
#' @aliases refFraction,ASEset-method
#' @docType methods
#' @param x \code{ASEset} object
#' @param ... arguments to forward to internal functions
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords reference fraction
#' @examples
#' 
#' #load example data
#' data(ASEset)
#' a <- ASEset
#' rf <- refFraction(a, strand="*")
#' 
#' @exportMethod refFraction
NULL

#' @rdname refFraction
setGeneric("refFraction", function(x, ... 
	){
    standardGeneric("refFraction")
})

setMethod("refFraction", signature(x = "ASEset"),
		function(x, strand="*",
			threshold.count.sample=0
	){
	
	#check for presence of genotype data
    if (!("genotype" %in% names(assays(x)))) {
		stop(paste("genotype matrix is not present as assay in",
				   " ASEset object, see '?inferGenotypes' "))
    }
	#check for presence of reference allele
	if(!("ref" %in% colnames(mcols(x)))){
		stop("column name 'ref' in mcols(x) is required")
	}

	#set all snp samples elements not heterozygote to zero
	acounts <- alleleCounts(x, strand=strand, return.class="array")
	acounts[array(!hetFilt(x), dim=c(nrow(x), ncol(x), 4))] <- 0

	#replace countmatrix in object x 
	x2 <-x
	alleleCounts(x2,strand=strand) <- acounts

	#calc frequency
	fr <- frequency(x2, strand=strand, return.class="array",
			threshold.count.sample = threshold.count.sample)

	#select only ref rows
	ar <- 	array(matrix(x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==mcols(x)[,"ref"]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x)))
	
	#rearrange to be able to transform back
	fr2 <- aperm(fr, c(3,2,1))
	ar2 <- aperm(ar, c(2,3,1))

	#rearrange dims
	refs <- fr2[ar2]
		
	matrix(refs,ncol=ncol(x),nrow=nrow(x), byrow=TRUE,dimnames=list(rownames(x),colnames(x)))


})


#' Reference allele
#' 
#' Extract the allele based on SNP location from the reference fasta file
#' 
#' The alleles will be placed in the rowData() meta column 'ref'
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
#' refAllele(ASEset.sim,fasta=fasta)
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
	
	ref <- scanFa(fl, param=rowData(x))
	close(fl)
	
	mcols(x)[["ref"]] <- as.vector(ref)
	x
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
setMethod("defaultMapBias", signature(x = "ASEset"), function(x){

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

		assays(x)[["mapBias"]] <- ar
		x
})



