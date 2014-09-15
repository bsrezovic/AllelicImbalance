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
#' @param strand strand option
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords reference fraction
#' @examples
#' 
#' #load example data
#' data(ASEset)
#' a <- ASEset
#'
#' #this example data contains to few SNPs to actually 
#' #measure the effect reliably, but as example it serves
#' #the purpose.
#'	
#' #prepare ASEset
#' genotype(a) <- inferGenotypes(a)
#' a <- refAllele(a,
#'		fasta=system.file('extdata/hg19.chr17.fa', 
#'		package='AllelicImbalance'))	
#'
#' rf <- refFraction(a, strand="*")
#' 
#' @exportMethod refFraction
NULL

#' @rdname refFraction
setGeneric("refFraction", function(x, strand="*" 
	){
    standardGeneric("refFraction")
})

setMethod("refFraction", signature(x = "ASEset"), function(x, strand="*"){
	
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
	alleleCounts(x,strand=strand) <- acounts

	#calc frequency
	fr <- frequency(x,strand=strand,return.class="array")

	#for loop in wait for a more ultimate solution
	ret <- matrix(NA, ncol=ncol(x), nrow=nrow(x),
			dimnames=list(rownames(x),colnames(x)))
	for (i in 1:nrow(x)){
		ret[i,] <- fr[i, , mcols(x)[i,"ref"] ]
	}
	
	ret
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
#' data(ASEset)
#'
#' fasta <- system.file('extdata/hg19.chr17.fa', package='AllelicImbalance')
#' refAllele(ASEset,fasta=fasta)
#' a <- refAllele(ASEset,fasta=fasta) 
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



