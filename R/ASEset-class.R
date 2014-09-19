#'@include initialize-methods.R
NULL

#' ASEset objects
#' 
#' Object that holds allele counts, genomic positions and map-bias for a set of
#' SNPs
#' 
#' An ASEset object differs from a regular SummarizedExperiment object in that
#' the assays contains an array instead of matrix. This array has ranges on the
#' rows, sampleNames on the columns and variants in the third dimension.
#' 
#' It is possible to use the commands barplot and locationplot on an ASEset
#' object see more details in \code{\link{barplot}} and
#' \code{\link{locationplot}}.
#' 
#' Three different alleleCount options are available. The simples one is the
#' * option, and is for experiments where the strand information is not
#' known e.g. non-stranded data. The unknown strand could also be for strand 
#' specific data when the aligner could not find
#' any strand associated with the read, but this should normally not happen, 
#' and if it does probably having an extremely low mapping quality. 
#' Then there are an option too add plus
#' and minus stranded data. When using this, it is essential to make sure that
#' the RNA-seq experiment under analysis has in fact been created so that
#' correct strand information was obtained. The most functions will by default 
#' have their strand argument set to '*'.
#' 
#' @name ASEset-class
#' @rdname ASEset-class
#' @aliases ASEset-class ASEset alleleCounts mapBias fraction arank 
#' frequency genotype genotype<- alleleCounts,ASEset-method mapBias,ASEset-method
#' fraction,ASEset-method arank,ASEset-method 
#' frequency,ASEset-method genotype,ASEset-method genotype<-,ASEset-method
#' alleleCounts<- alleleCounts<-,ASEset-method
#' 
#' @docType class
#' @param x ASEset object
#' @param strand which strand of '+', '-' or '*'
#' @param verbose makes function more talkative
#' @param return.type return 'names', rank or 'counts'
#' @param return.class return 'list' or 'array'
#' @param value value as replacement
#' @param ... additional arguments
#' @return An object of class ASEset containing location information and allele
#' counts for a number of SNPs measured in a number of samples on various
#' strand, as well as mapBias information. All data is stored in a manner
#' similar to the \code{\link[GenomicRanges]{SummarizedExperiment}} class.
#' @section Table: table(...)
#' 
#' \describe{
#' Arguments: \item{...}{An \code{ASEset object} that contains the
#' variants of interest} 
#'
#' The generics for table does not easily allow more than one argument
#' so in respect to the different strand options, \code{table} will
#' return a SimpleList with length 3, one element for each strand.
#' }
#'
#' @section Frequency: frequency(x, 
#' return.class = "list", strand = "*",
#' threshold.count.sample = 15)
#' 
#' \describe{
#' Arguments: \item{x}{An \code{ASEset object} that contains the
#' variants of interest} 
#' 
#' \item{x}{threshold.count.samples} if sample has fewer counts the function
#' return NA.
#'	
#' }
#' @section Constructor: ASEsetFromCountList(rowData, countListNonStranded =
#' NULL, countListPlus = NULL, countListMinus = NULL, countListUnknown = NULL,
#' colData = NULL, mapBiasExpMean = array(),verbose=FALSE ...)
#' 
#' \describe{
#' 
#' Arguments: \item{rowData}{A \code{GenomicRanges object} that contains the
#' variants of interest} \item{countListNonStranded}{A \code{list} where each
#' entry is a matrix with allele counts as columns and sample counts as rows}
#' \item{countListPlus}{A \code{list} where each entry is a matrix with allele
#' counts as columns and sample counts as rows} \item{countListMinus}{A
#' \code{list} where each entry is a matrix with allele counts as columns and
#' sample counts as rows} \item{countListUnknown}{A \code{list} where each
#' entry is a matrix with allele counts as columns and sample counts as rows}
#' \item{colData}{A \code{DataFrame} object containing sample specific data}
#' \item{mapBiasExpMean}{A 3D \code{array} describing mapping bias. The SNPs
#' are in the 1st dimension, samples in the 2nd dimension and variants in the
#' 3rd dimension.} \item{verbose}{Makes function more talkative}
#' \item{...}{arguments passed on to SummarizedExperiment constructor} }
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The
#' \code{\link[GenomicRanges]{SummarizedExperiment}} for ranges operations.  }
#' @keywords class ASEset
#' @examples
#' 
#' 
#' #make example countList
#' set.seed(42)
#' countListPlus <- list()
#' snps <- c('snp1','snp2','snp3','snp4','snp5')
#' for(snp in snps){
#'   count<-matrix(rep(0,16),ncol=4,dimnames=list(
#' c('sample1','sample2','sample3','sample4'),
#' c('A','T','G','C')))
#'   #insert random counts in two of the alleles 
#'   for(allele in sample(c('A','T','G','C'),2)){
#' count[,allele]<-as.integer(rnorm(4,mean=50,sd=10))
#'   }
#'   countListPlus[[snp]] <- count
#' }
#' 
#' #make example rowData
#' rowData <- GRanges(
#'   seqnames = Rle(c('chr1', 'chr2', 'chr1', 'chr3', 'chr1')),
#'   ranges = IRanges(1:5, width = 1, names = head(letters,5)),
#'   snp = paste('snp',1:5,sep='')
#' )
#'
#' #make example colData
#' colData <- DataFrame(Treatment=c('ChIP', 'Input','Input','ChIP'), 
#'  row.names=c('ind1','ind2','ind3','ind4'))
#' 
#' #make ASEset 
#' a <- ASEsetFromCountList(rowData, countListPlus=countListPlus, 
#' colData=colData)
#' 
#'
#' @exportClass ASEset
#' @exportMethod alleleCounts alleleCounts<- mapBias fraction arank
#' frequency genotype genotype<-
#' 
#' @export 

setClass("ASEset", contains = "SummarizedExperiment", 
	representation(variants = "vector"))

#' @rdname ASEset-class
setGeneric("alleleCounts", function(x, strand = "*", return.class="list") {
    standardGeneric("alleleCounts")
})

setMethod("alleleCounts", signature(x = "ASEset"), function(x, strand = "*",
	return.class="list") {

    if (!sum(strand %in% c("+", "-", "*")) > 0) {
        stop("strand parameter has to be either '+', '-', '*' ")
    }
    
    if (strand == "+") {
        el <- "countsPlus"
    } else if (strand == "-") {
        el <- "countsMinus"
    } else if (strand == "*") {
        el <- "countsUnknown"
    } else {
        stop("not existing strand option")
    }
    
    # check if strand option is present as assay
    if (!(el %in% names(assays(x)))) {
		stop(paste("neither '+' '-' or '*'strand is present as assay in",
				   " ASEset object"))
    }
    
	if(return.class=="array"){
		if(el=="combine"){
			ar <- assays(x)[["countsPlus"]] + assays(x)[["countsMinus"]]
		}else{
			ar <- assays(x)[[el]]
		}
		#add variant names (preferably this could take place during 
		# initialization of the ASEset object)
		dimnames(ar)[[3]]<- x@variants
		ar
	}else if(return.class=="list"){

		if(el=="combine"){
			ar <- assays(x)[["countsPlus"]] + assays(x)[["countsMinus"]]
		}else{
			ar <- assays(x)[[el]]
		}

		alleleCountList <- list()

		for (i in 1:nrow(ar)) {
			mat <- ar[i, , ]

			if (class(mat) == "integer") {
				mat <- t(as.matrix(mat))
			}
			if (class(mat) == "numeric") {
				mat <- t(mat)
			}
			colnames(mat) <- x@variants
			rownames(mat) <- colnames(x)

			alleleCountList[[i]] <- mat
		}

		# add snp id
		names(alleleCountList) <- rownames(x)
		
		# return object
		alleleCountList
		
	}else{
		stop("return.class has to be 'list' or 'array'")
	}

    
})

#' @rdname ASEset-class
setGeneric("alleleCounts<-", function(x, strand = "*", value) {
    standardGeneric("alleleCounts<-")
})

setMethod("alleleCounts<-", signature(x = "ASEset"), function(x,
	 strand = "*", value) {

    if (!sum(strand %in% c("+", "-", "*")) > 0) {
        stop("strand parameter has to be either '+', '-', '*' ")
    }
    
    if (strand == "+") {
        el <- "countsPlus"
    } else if (strand == "-") {
        el <- "countsMinus"
    } else if (strand == "*") {
        el <- "countsUnknown"
    } else {
        stop("not existing strand option")
    }
    
	#check that value has the right dimensions
	if(!all(dim(value)==dim(assays(x)[["countsUnknown"]]))){
		stop("dimensions in replacement object is not correct")
	}

	assays(x)[[el]] <- value
	x

})

#' @rdname ASEset-class
setGeneric("mapBias", function(x, ...) {
    standardGeneric("mapBias")
})

setMethod("mapBias", signature(x = "ASEset"), function(x,
	return.class="list") {
    # assume alleleCount information is stored as element 1

	if(return.class=="array"){
		return(assays(x)[["mapBias"]])

	}else if(return.class=="list"){
		mapBiasList <- list()
		for (i in 1:nrow(x)) {
			mat <- assays(x)[["mapBias"]][i, , ]
			if (class(mat) == "numeric") {
				dim(mat) <- c(1, 4)
				rownames(mat) <- colnames(x)
			}
			colnames(mat) <- x@variants
			
			mapBiasList[[i]] <- mat
		}
		# add snp id
		names(mapBiasList) <- rownames(x)
		
		return(mapBiasList)
	}
    
})

#' @rdname ASEset-class
setGeneric("fraction", function(x, strand = "*", verbose = FALSE) {
    standardGeneric("fraction")
})

setMethod("fraction", signature(x = "ASEset"), function(x, strand = "*", 
    verbose = FALSE) {
    
    if (!sum(strand %in% c("+", "-", "*")) > 0) {
        stop("strand parameter has to be either '+', '-', '*' ")
    }
    
    if (strand == "+") {
        el <- "countsPlus"
    } else if (strand == "-") {
        el <- "countsMinus"
    } else if (strand == "*") {
        el <- "countsUnknown"
    } else {
        stop("unknown strand option")
    }
    
    # check if strand option is present as assay
    if (!(el %in% names(assays(x)))) {
        stop(paste("strand", strand, " is not present as assay in ASEset object"))
    }
    
    fractionList <- list()
    
    for (i in 1:nrow(x)) {
        # getting revelant data
        tmp <- alleleCounts(x, strand)[[i]]
        
        # calculating major and minor allele, warn if the two remaining alleles have too
        # many counts
        if (nrow(tmp) > 1) {
            countsByAllele <- apply(tmp, 2, sum, na.rm = TRUE)
        } else {
            countsByAllele <- tmp
        }
        majorAllele <- colnames(tmp)[order(countsByAllele, decreasing = TRUE)][1]
        minorAllele <- colnames(tmp)[order(countsByAllele, decreasing = TRUE)][2]
        majorAndMinorFraction <- sum(countsByAllele[c(majorAllele, minorAllele)])/sum(countsByAllele)
        if (verbose & majorAndMinorFraction < 0.9) {
            cat(paste("Snp", "was possible tri-allelic, but only two most",
					  "frequent alleles were plotted. Counts:"), 
					"\n")
            cat(paste(paste(names(countsByAllele), countsByAllele, sep = "="), collapse = ", "), 
                "\n")
        }
        
        # calculating percentage and ylim and setting no-count samples to colour grey90
        fraction <- tmp[, majorAllele]/(tmp[, majorAllele] + tmp[, minorAllele])
        # fraction[is.nan(fraction)]<-1
        fractionList[[i]] <- fraction
        
        if (i%%300 == 0) {
            cat(paste("processed ", i, " snps\n", sep = ""))
        }
    }
    
    names(fractionList) <- rownames(x)
    
    # return object fractionList
    m <- as.matrix(as.data.frame(fractionList))
    rownames(m) <- colnames(x)
    
    m
    
})

#' @rdname ASEset-class
setGeneric("arank", function(x, return.type = "names", 
	return.class = "list", strand = "*", ...) {
    standardGeneric("arank")
})

setMethod("arank", signature(x = "ASEset"), function(x, return.type = "names", 
	return.class = "list", strand = "*", ...) {
	
	if(return.class=="matrix"){
		if (return.type == "names") {
			ar <- t(apply(apply(alleleCounts(x,
						strand=strand,return.class="array"),c(1,3),sum),
					1, function(x){rank(x,ties.method="first")}))
			ar2 <- t(apply(ar,1,function(x){
						   x <- sort(x,index.return=TRUE,decreasing=TRUE)$ix
						   x
					}))

			mat <- matrix(x@variants[ar2],ncol=4, nrow=nrow(x), byrow=FALSE,
					  dimnames=list(dimnames(ar)[[1]],c(1,2,3,4)))
		
			return(mat)
		}else if (return.type == "rank") {
			return(t(apply(apply(alleleCounts(x,
						strand=strand,return.class="array"),c(1,3),sum),
					1, function(x){rank(x,ties.method="first")})))
		}else if (return.type == "counts") {
			return(apply(alleleCounts(x,strand=strand,return.class="array"),c(1,3),sum))
		}else{stop("return.type is not valid")}

	}

	if(return.class=="list"){
		acounts <- alleleCounts(x, strand = strand)
		
		if (return.type == "names") {
			lapply(acounts, function(x) {
				names(sort(apply(x, 2, sum), decreasing = TRUE))
			})
		}else if (return.type == "counts") {
			lapply(acounts, function(x) {
				as.numeric(sort(apply(x, 2, sum), decreasing = TRUE))
			})
		}else if (return.type == "rank") {
			stop("rank for return.class 'list' is not implemented yet")
		}else{stop("return.type is not valid")}
	}
}) 

#setGeneric("table")
#setGeneric("table", function(x, strand = "*", sortBy="none", ...) {
#    standardGeneric("table")
#})

# @rdname ASEset-class
#setMethod("table", signature(... = "ASEset"), function(...) {
#
#	args <- list(...)
#	if (length(args) > 1)
#	  stop("Only one argument in '...' supported")
#	x <- args[[1L]]
#
#	#because the generis of table is rubbish we have to return a list for each strand
#	#retList <- list()
#
#	#for(strand in c("+","-","*")){
#	#	df <- data.frame(row.names=rownames(x))
#	#	df[,c("chromosome","position")] <- c(as.character(seqnames(x)),start(x))
#	#	df <- cbind(df,as.data.frame(arank(x, return.type="counts",
#	#				   return.class="matrix",strand=strand)))
#
#	#	#if only one sample add fraction to table()
#	#	if(ncol(x)==1){
#	#		df[,"fraction"] <- as.vector(fraction(x,strand=strand))
#	#	}
#
#	#	retList[[strand]] <- df
#
#	}
#	#return(SimpleList(retList))
#
#	df()
#
#})	

#' @rdname ASEset-class
setGeneric("frequency")

setMethod("frequency", signature(x = "ASEset"), function(x, 
	return.class = "list", strand = "*",
	threshold.count.sample = 15) {

	ar <- alleleCounts(x, strand=strand, return.class="array")
	allele.count.tot <- apply(ar, c(1,2), sum)

	#set allele.count.tot to NaN for all values not passing threshold
	#this is to indicate that we do not have enough reads to say anything
	tf <- allele.count.tot  < threshold.count.sample
	allele.count.tot[tf] <- NaN

	#make frequency array
	ar <- ar * array(as.vector(1/allele.count.tot),dim=dim(ar))

	if(return.class=="array"){
		return(ar)
	}else if(return.class=="list"){
		lst <- list()
		for (i in 1:nrow(x)){
			mat <- ar[i,,]
			dimnames(mat) <- list(colnames(x),x@variants)
			lst[[i]] <- mat
		}
		names(lst) <- rownames(x)
		lst
	}else{
		stop("return.class has to be 'array' or 'list'")
	}
})

#' @rdname ASEset-class
setGeneric("genotype", function(x){
    standardGeneric("genotype")
})

setMethod("genotype", signature(x = "ASEset"), function(x){
	
    if (!("genotype" %in% names(assays(x)))) {
		stop(paste("genotype matrix is not present as assay in",
				   " ASEset object, see '?inferGenotypes' "))
    }
	assays(x)[["genotype"]]
})
#' @rdname ASEset-class
setGeneric("genotype<-", function(x,value){
    standardGeneric("genotype<-")
})

setMethod("genotype<-", signature(x = "ASEset"), function(x,value){
	
	#check dimensions
	if(!nrow(x)==nrow(value)){
		stop("nrow(x) is not equal to nrow(value)")	
	}
	if(!ncol(x)==ncol(value)){
		stop("ncol(x) is not equal to ncol(value)")	
	}

	assays(x)[["genotype"]] <- value	
	x
})

#' @rdname ASEset-class
setGeneric("countsPerSnp", function(x, ...){
    standardGeneric("countsPerSnp")
})

#could be renamed to countsAllAlleles
setMethod("countsPerSnp", signature(x = "ASEset"), function(x, 
	return.class = "matrix", return.type="mean", strand = "*") {

	if(return.class=="matrix"){
		return(apply(alleleCounts(x, strand=strand, return.class="array"), c(1,2), sum))
	}else if(return.class=="vector"){
		if(return.type=="all"){
			return(apply(alleleCounts(x, strand=strand, return.class="array"), 1, sum))
		}else if(return.type=="mean"){
			return(apply(apply(alleleCounts(x, strand=strand, return.class="array"), c(1,2), sum), 1, mean))
		}
	}else{
		stop("return.class has to be 'vec'")
	}
})

