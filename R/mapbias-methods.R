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
			threshold.count.sample=1,
			random.ref=FALSE,
			...
	){
	
	#check for presence of genotype data
    if (!("genotype" %in% names(assays(x)))) {
		stop(paste("genotype matrix is not present as assay in",
				   " ASEset object, see '?inferGenotypes' "))
    }
	#check for presence of reference allele
	if(!random.ref){
		if(!("ref" %in% colnames(mcols(x)))){
			stop("column name 'ref' in mcols(x) is required")
		}
	}else{
		mcols(x)[,"ref"] <- randomRef(x, inferGenotypes=TRUE)
	}

	#set all snp samples elements not heterozygote to zero
	acounts <- alleleCounts(x, strand=strand, return.class="array")
	acounts[array(!hetFilt(x), dim=c(nrow(x), ncol(x), 4))] <- 0

	#replace countmatrix in object x 
	x2 <-x
	alleleCounts(x2,strand=strand) <- acounts

	#calc frequency threshold.count.sample vector will be used further down.
	#not using it here is because we dont want to calculate this many times.
	#because it is heavy and not needed to calculate more than one time.
	fr <- frequency(x2, strand=strand, return.class="array",
			threshold.count.sample = 1)

	#select only ref rows
	ar <- 	array(matrix(x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==mcols(x)[,"ref"]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x)))
	
	#rearrange to be able to transform back
	fr2 <- aperm(fr, c(3,2,1))
	ar2 <- aperm(ar, c(2,3,1))


	#subset ref allele frequencies
	refs <- fr2[ar2]

	#prepare allele.count.tot to test against the threshold.count.sample
	allele.count.tot <- apply(acounts, c(1,2), sum)

	ar.tot  <- array(allele.count.tot,dim=c(nrow(x),ncol(x),length(threshold.count.sample))) < 
					aperm(array(threshold.count.sample,
						dim=c(length(threshold.count.sample),nrow(x),ncol(x) )), 
						c(2,3,1))
	ar.ref <- aperm(array(refs, dim=c(ncol(x),nrow(x),length(threshold.count.sample)),
						  dimnames=list(colnames(x),rownames(x),paste(threshold.count.sample))),
					c(2,1,3))

	#set hte ones not passing threshold to NA
	ar.ref[ar.tot] <- NA

	#}
		
	#return object
	ar.ref

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

#' Random ref allele from genotype
#' 
#' Create a vector of random reference alleles
#' 
#' Randomly shuffles which of the two alleles for each genotype that is 
#' indicated as reference allele.
#' 
#' @name randomRef
#' @rdname randomRef
#' @aliases randomRef,ASEset-method
#' @docType methods
#' @param x \code{ASEset} object
#' @param ... internal arguments
#' @param inferGenotypes infer genotypes from count matrix
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords mapbias
#' @examples
#' 
#' #load example data
#' data(ASEset.sim)
#' a <- ASEset.sim
#'
#' mcols(a)[["ref"]] <- randomRef(a, inferGenotypes=TRUE) 
#'
NULL

#' @rdname randomRef
#' @export
setGeneric("randomRef", function(x,... ){
    standardGeneric("randomRef")
})

#' @rdname randomRef
#' @export
setMethod("randomRef", signature(x = "ASEset"), function(x, inferGenotypes=FALSE){

		if (!inferGenotypes){
			if (!("genotype" %in% names(assays(x)))) {
				stop(paste("genotype matrix is not present as assay in",
						   " ASEset object, see '?inferGenotypes' "))
			}
		}else{
			genotype(x) <- inferGenotypes(x)
		}

		g <- genotype(x)
		allele1 <- sub("/.*","",as.character(g))
		allele2 <- sub(".*/","",as.character(g))

		mat1 <- matrix(allele1, nrow=nrow(x),ncol=ncol(x),byrow=FALSE)
		mat2 <- matrix(allele2, nrow=nrow(x),ncol=ncol(x),byrow=FALSE)

		mat <- matrix(c(mat1,mat2),nrow=nrow(x))

		alleles <- apply(mat, 1, 
					 function(x){names(sort(table(x), decreasing=TRUE))[c(1,2)]})

		apply(alleles, 2, sample, size=1)

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


