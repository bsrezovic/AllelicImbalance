#'@include AllelicImbalance-package.R
NULL

# setMethods('initialize') is not required at the moment

#' Initialize ASEset
#' 
#' Functions to construct ASEset objects
#' 
#' The resulting ASEset object is based on the RangedSummarizedExperiment
#' class, and will therefore inherit the same accessors and ranges operations.
#' 
#' If both countListPlus and countListMinus are given they will be used to 
#' calculate countListUnknown, which is the sum of the plus and minus strands.
#' 
#' countListPlus, countListMinus and countListUnknown are
#' i.e. the outputs from the getAlleleCounts function.
#'
#' aquals is new for the devel branch and will be changed slighly before the relase 
#' to include better granularity.
#' 
#' @name initialize-ASEset
#' @rdname initialize-ASEset
#' @aliases initialize-ASEset ASEsetFromCountList
#' @param rowRanges A \code{GenomicRanges object} that contains the variants of
#' interest
#' @param countListPlus A \code{list} where each entry is a matrix with allele
#' counts as columns and sample counts as rows
#' @param countListMinus A \code{list} where each entry is a matrix with allele
#' counts as columns and sample counts as rows
#' @param countListUnknown A \code{list} where each entry is a matrix with
#' allele counts as columns and sample counts as rows
#' @param countsPlus An array containing the countinformation
#' @param countsMinus An array containing the countinformation
#' @param countsUnknown An array containing the countinformation
#' @param aquals A 4-D array containing the countinformation, see details
#' @param genotype matrix
#' @param colData A \code{DataFrame} object containing sample specific data
#' @param phase A \code{matrix} or an \code{array} containing phase information. 
#' @param mapBiasExpMean A 3D \code{array} where the SNPs are in the 1st
#' dimension, samples in the 2nd dimension and variants in the 3rd dimension.
#' @param verbose Makes function more talkative
#' @param ... arguments passed on to SummarizedExperiment constructor
#' @return \code{ASEsetFromCountList} returns an \code{ASEset} object.
#' @note \code{ASEsetFromCountList} requires the same input data as a
#' RangedSummarizedExperiment, but with minimum one assay for the allele counts.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item
#' \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} objects. }
#' @keywords ASEset ASEsetFromCountList
#' @examples
#' 
#' #make example alleleCountListPlus
#' set.seed(42)
#' countListPlus <- list()
#' snps <- c('snp1','snp2','snp3','snp4','snp5')
#' for(snp in snps){
#' count<-matrix(rep(0,16),ncol=4,dimnames=list(
#' c('sample1','sample2','sample3','sample4'),
#' c('A','T','G','C')))
#' #insert random counts in two of the alleles 
#' for(allele in sample(c('A','T','G','C'),2)){
#' count[,allele]<-as.integer(rnorm(4,mean=50,sd=10))
#' }
#' countListPlus[[snp]] <- count
#' }
#' 
#' #make example alleleCountListMinus
#' countListMinus <- list()
#' snps <- c('snp1','snp2','snp3','snp4','snp5')
#' for(snp in snps){
#' count<-matrix(rep(0,16),ncol=4,dimnames=list(
#' c('sample1','sample2','sample3','sample4'),
#' c('A','T','G','C')))
#' #insert random counts in two of the alleles 
#' for(allele in sample(c('A','T','G','C'),2)){
#' count[,allele]<-as.integer(rnorm(4,mean=50,sd=10))
#' }
#' countListMinus[[snp]] <- count
#' }
#' 
#' 
#' #make example rowRanges
#' rowRanges <- GRanges(
#' seqnames = Rle(c('chr1', 'chr2', 'chr1', 'chr3', 'chr1')),
#'          ranges = IRanges(1:5, width = 1, names = head(letters,5)),
#'          snp = paste('snp',1:5,sep='')
#'          )
#' #make example colData
#' colData <- DataFrame(Treatment=c('ChIP', 'Input','Input','ChIP'), 
#'  row.names=c('ind1','ind2','ind3','ind4'))
#' 
#' #make ASEset 
#' a <- ASEsetFromCountList(rowRanges, countListPlus=countListPlus, 
#' countListMinus=countListMinus, colData=colData)
#' 
NULL


#' @rdname initialize-ASEset
#' @export 
ASEsetFromCountList <- function(rowRanges, countListUnknown = NULL, countListPlus = NULL, 
    countListMinus = NULL, colData = NULL, mapBiasExpMean = NULL, phase = NULL, aquals = NULL,
    verbose = FALSE, ...) {
    
    if (verbose) {
        cat("rowRanges\n")
        cat(class(rowRanges))
        cat("countListPlus\n")
        cat(class(countListPlus))
        cat("countListMinus\n")
        cat(class(countListMinus))
        cat("countListUnknown\n")
        cat(class(countListUnknown))
        cat("colData\n")
        cat(class(colData))
        cat("mapBiasExpMean\n")
        cat(class(mapBiasExpMean))
    }

    # check that at least one of the countList options are not null
    if (is.null(c(countListPlus, countListMinus, countListUnknown))) {
        stop("at least one of the countList options has to be specified")
    }
    
    countLists <- c("countListPlus", "countListMinus", "countListUnknown")[(c( 
        !is.null(countListPlus), !is.null(countListMinus), !is.null(countListUnknown)))]
    
    # check that all lengths are the same in all lists
    l <- unlist(lapply(countLists, function(x) {
        length(get(x))
    }))
    if (sum(!(l == l[1])) > 0) {
        stop("one or more list contains more or less list elements than the others")
    }
    
    # check that the number of columns in each dataframe are the same
    l <- unlist(lapply(countLists, function(x) {
        lapply(get(x), ncol)
    }))
    mat <- matrix(unlist(l), ncol = length(l))
    if (sum(!(mat == mat[1])) > 0) {
        stop("one or more list contains more or less columns than the others")
    }
    
    # check that the number of rows in each dataframe are the same
    l <- unlist(lapply(countLists, function(x) {
        lapply(get(x), nrow)
    }))
    mat <- matrix(unlist(l), ncol = length(l))
    if (sum(!(mat == mat[1])) > 0) {
        stop("one or more list contains more or less rows than the others")
    }
    
    # check that all colnames in all dataframes are the same in same order
    cMatrix <- matrix(NA, nrow = length(unlist(unique(lapply(countLists[1], colnames)))), 
        ncol = length(countLists))
    for (i in 1:length(countLists)) {
        # within list comparision
        countListName <- countLists[i]
        countList <- get(countListName)
        l <- lapply(countList, colnames)
        mat2 <- matrix(unlist(l), ncol = length(l))
        if (sum(!(mat2 == mat2[, 1])) > 0) {
            stop(paste("the names or order of names are not the same in all data frames in list", 
                countListName))
        }
        cMatrix[, i] <- mat2[, 1]
    }
    if (sum(!(cMatrix == cMatrix[, 1])) > 0) {
        stop(paste("the names or order of names are not the same between all data frame lists"))
    }
    
    # check that all rownames in all dataframes are the same in same order
    cMatrix <- matrix(NA, nrow = length(unlist(unique(lapply(countLists[1], rownames)))), 
        ncol = length(countLists))
    for (i in 1:length(countLists)) {
        # within list comparision
        countListName <- countLists[i]
        countList <- get(countListName)
        l <- lapply(countList, rownames)
        mat2 <- matrix(unlist(l), ncol = length(l))
        if (sum(!(mat2 == mat2[, 1])) > 0) {
            stop(paste("the names or order of names are not the same in all data frames in list", 
                countListName))
        }
        cMatrix[, i] <- mat2[, 1]
    }
    if (sum(!(cMatrix == cMatrix[, 1])) > 0) {
        stop(paste("the names or order of names are not the same between all data frame lists"))
    }
    
    # check mapBiasExpMean
    if (!(is.null(mapBiasExpMean))) {
        if (!class(mapBiasExpMean) == "array") {
            stop("mapBiasExpMean has to be of class array")
        }
        if (!length(dim(mapBiasExpMean)) == 3) {
            stop("mapBiasExpMean has to have three dimensions")
        }
	}
        
    # choose a common countList by picking the first one, for dimension info
    countList <- get(countLists[1])
    ind <- length(unlist(unique(lapply(countList, rownames))))
    snps <- length(countList)
    
    # SimpleList init
    assays <- SimpleList()
    
    # plus
    if (!is.null(countListPlus)) {
		#empty array that handles only four nucleotides 
        ar1 <- array(NA, c(snps, ind, 4))  
        for (i in 1:snps) {
            ar1[i, , ] <- countListPlus[[i]]
        }
    }
    # minus
    if (!is.null(countListMinus)) {
		#empty array that handles only four nucleotides 
        ar2 <- array(NA, c(snps, ind, 4))  
        for (i in 1:snps) {
            ar2[i, , ] <- countListMinus[[i]]
        }
        
    }
    # unknown
    if (!is.null(countListUnknown)) {
		#empty array that handles only four nucleotides 
        ar3 <- array(NA, c(snps, ind, 4))  
        for (i in 1:snps) {
            ar3[i, , ] <- countListUnknown[[i]]
        }

	}	

	#ar <- array(0, c(snps, ind, 4, 2))  
	if((is.null(countListMinus)) & (is.null(countListPlus))){
	#	ar[,,,1] <- ar3
		assays[["countsPlus"]] <- ar3
	}else{
		if(!is.null(countListPlus)){
	#		ar[,,,1] <- ar1	
			assays[["countsPlus"]] <- ar1
		}
		if(!is.null(countListMinus)){
	#		ar[,,,2] <- ar2	
			assays[["countsMinus"]] <- ar2
		}
	}
	#assays[["acounts"]] <- ar

    # assign mapBiasExpMean
    if (is.null(mapBiasExpMean)) {
        assays[["mapBias"]] <- getDefaultMapBiasExpMean3D(countList)
    } else {
        assays[["mapBias"]] <- mapBiasExpMean
    }
    
    # assign phase if user provides it
    if (is.null(phase)) {
        assays[["phase"]] <- defaultPhase(snps,ind)
    }

    if (is.null(aquals)) {
        assays[["aquals"]] <- NULL
    }

    if (is.null(colData)) {
        colData <- DataFrame(row.names = unlist(unique(lapply(countList, rownames))))
    }
    
    sset <- SummarizedExperiment(assays = assays, rowRanges = rowRanges, colData = colData, 
        ...)
    
    rownames(sset) <- names(countList)
    
    
    # use colnames in list matrices as variants
    variants <- unlist(unique(lapply(countList, colnames)))
    
    ASEset <- function(sset, variants) {
        # create object
        new("ASEset", sset, variants = variants)
    }
    
    # create object
    ASEset(sset, variants = variants)
} 

#' @rdname initialize-ASEset
#' @export 
ASEsetFromArrays <- function(rowRanges, countsUnknown = NULL, countsPlus = NULL, 
    countsMinus = NULL, colData=NULL,  mapBiasExpMean = NULL, phase = NULL,
	genotype = NULL, aquals = NULL,
    verbose = FALSE, ...) {
    
    # check that at least one of the countList options are not null
    if (is.null(c(countsPlus, countsMinus, countsUnknown))) {
        stop("at least one of the countList options has to be specified")
    }
    
    countLists <- c("countsPlus", "countsMinus", "countsUnknown")[(c( 
        !is.null(countsPlus), !is.null(countsMinus), !is.null(countsUnknown)))]
    
    # check mapBiasExpMean
    if (!(is.null(mapBiasExpMean))) {
        if (!class(mapBiasExpMean) == "array") {
            stop("mapBiasExpMean has to be of class array")
        }
        if (!length(dim(mapBiasExpMean)) == 3) {
            stop("mapBiasExpMean has to have three dimensions")
        }
	}    

    # choose a common countList by picking the first one, for dimension info
    countList <- get(countLists[1])
    ind <- dim(countList)[2]
    ind.names <- dimnames(countList)[2]
    snps <- dim(countList)[1]
    snps.names <- dimnames(countList)[1]
    
    # SimpleList init
    assays <- SimpleList()
    
	ar <- array(0, c(snps, ind, 4, 2)) 

	if((is.null(countsMinus)) & (is.null(countsPlus))){
		#ar[,,,1] <- countsUnknown
		assays[["countsPlus"]] <- countsUnknown
	}else{
		if(!is.null(countsPlus)){
			#ar[,,,1] <- countsPlus	
		assays[["countsPlus"]] <- countsPlus
		}
		if(!is.null(countsMinus)){
			#ar[,,,2] <- countsMinus	
		assays[["countsMinus"]] <- countsMinus
		}
	}

    # assign mapBiasExpMean
    if (is.null(mapBiasExpMean)) {
        assays[["mapBias"]] <- getDefaultMapBiasExpMean3D(ar)
    } else {
        assays[["mapBias"]] <- mapBiasExpMean
    }
    
    # assign phase if user provides it
    if (is.null(phase)) {
        assays[["phase"]] <- defaultPhase(snps,ind)
    }else {
        assays[["phase"]] <- phase
	}
    if (!is.null(genotype)) {
        assays[["genotype"]] <- genotype
    }

    if (is.null(aquals)) {
        assays[["aquals"]] <- NULL
    }

    if (is.null(colData)) {
		sset <- SummarizedExperiment(assays = assays, rowRanges = rowRanges,
        ...)
    }else{
		sset <- SummarizedExperiment(assays = assays, rowRanges = rowRanges, colData=colData,
        ...)
	}

    
    rownames(sset) <- rownames(countList)
    
    
    # use colnames in list matrices as variants
    #variants <- dimnames(assays[["countsUnknown"]])[[3]]    
    variants <- c( "A",  "C", "G", "T")
    ASEset <- function(sset, variants) {
        # create object
        new("ASEset", sset, variants = variants)
    }
    
    # create object
    ASEset(sset, variants = variants)
} 

#' Initialize ReferenceBias
#' 
#' Functions to construct ReferenceBias objects
#' 
#' produces a class container for reference bias calculations
#' 
#' @name initialize-ReferenceBias
#' @rdname initialize-ReferenceBias
#' @aliases initialize-ReferenceBias RBias
#' @param x \code{ASEset} 
#' @param ... internal arguments
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords bias mapbias refBias
#' @examples
#'
#' data(ASEset)
#' a <- ASEset
#' refbiasObject <- RBias(a)
#' 
NULL


#' @rdname initialize-ReferenceBias
#' @export 
#setMethod("ReferenceBias","ReferenceBias", function(
RBias <- function(
	x = "ASEset",
	...
	){
	#if non-stranded data	
	if(all(c("countsPlus","countsMinus") %in% names(assays(x)))){
		assay <- 
			array(c(refFraction(x,strand="*",...),
				  refFraction(x,strand="+",...),
				  refFraction(x,strand="-",...)),
				  dim=c(nrow(x),ncol(x),3),
				  dimnames=list(rownames(x),colnames(x),c("*","+","-")))
	}
	else if(c("countsUnknown") %in% names(assays(x))){
		assay <- 
			array(c(refFraction(x,strand="*",...),
				  matrix(NA, nrow=nrow(x),ncol=ncol(x)),
				  matrix(NA, nrow=nrow(x),ncol=ncol(x))),
				  dim=c(nrow(x),ncol(x),3),
				  dimnames=list(rownames(x),colnames(x),c("*","+","-")))
	}

	sset <- SummarizedExperiment(assays = SimpleList(referenceFrequency=assay), rowRanges = rowRanges(x), colData = colData(x)) 
	rownames(sset) <- rownames(x)

	#valid
	#validObject(.Object)

	#Return object
	new("ReferenceBias", sset, strands = c("*","+","-"))
}

#' Initialize DetectedAI
#' 
#' Functions to construct DetectedAI objects
#' 
#' produces a class container for reference bias calculations
#' 
#' @name initialize-DetectedAI
#' @rdname initialize-DetectedAI
#' @aliases initialize-DetectedAI 
#' @param x \code{ASEset} 
#' @param strand set strand to detectAI over "+","-","*"
#' @param reference.frequency frequencies of reference alleles based allele counts
#' @param threshold.frequency logical array for frequency thresholds
#' @param threshold.count.sample logical array for per sample allele count thresholds
#' @param threshold.delta.frequency logical array for delta frequency thresholds.
#' @param threshold.pvalue logical array for pvalue thresholds (max 1, min 0)
#' @param threshold.frequency.names character vector
#' @param threshold.count.sample.names character vector
#' @param threshold.delta.frequency.names character vector
#' @param threshold.pvalue.names character vector
#' @param ... internal arguments
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords bias mapbias refBias
#' @examples
#'
#' data(ASEset)
#' a <- ASEset
#' dai <- detectAI(a)
#' 
NULL

#' @rdname initialize-DetectedAI
#' @export 
#setMethod("DetectedAI","DetectedAI", function(
DetectedAIFromArray <- function(
	x = "ASEset",
	strand = "*",
	reference.frequency=NULL,
	threshold.frequency=NULL,
	threshold.count.sample=NULL,
	threshold.delta.frequency=NULL,
	threshold.pvalue=NULL,

	#reference.frequency.names=NULL,
	threshold.frequency.names=NULL,
	threshold.count.sample.names=NULL,
	threshold.delta.frequency.names=NULL,
	threshold.pvalue.names=NULL,

	...){

	sset <- SummarizedExperiment(
				assays = SimpleList(
					reference.frequency=reference.frequency,
					threshold.frequency=threshold.frequency,
					threshold.count.sample=threshold.count.sample,
					threshold.delta.frequency=threshold.delta.frequency,
					threshold.pvalue=threshold.pvalue
					), 
				rowRanges = rowRanges(x), 
				colData = colData(x)
			) 

	rownames(sset) <- rownames(x)

	#valid
	#validObject(.Object)

	#Return object
	new("DetectedAI", sset,
		strand = strand,
		#reference.frequency.names=reference.frequency.names,
		threshold.frequency.names=threshold.frequency.names,
		threshold.count.sample.names=threshold.count.sample.names,
		threshold.delta.frequency.names=threshold.delta.frequency.names,
		threshold.pvalue.names=threshold.pvalue.names
	)
}

###################
#
# Global reference frequency analysis
#
###################

#' Initialize GlobalAnalysis
#' 
#' Functions to construct GlobalAnalysis objects
#' 
#' produces a class container for a global analysis
#' 
#' @name initialize-GlobalAnalysis
#' @rdname initialize-GlobalAnalysis
#' @aliases initialize-GlobalAnalysis
#' @param x \code{ASEset} 
#' @param ... internal arguments
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords global
#' @examples
#'
#' data(ASEset)
#' a <- ASEset
#' # gba <- gba(a)
#' 
NULL


#' @rdname initialize-GlobalAnalysis
#' @export 
#setMethod("ReferenceBias","ReferenceBias", function(
GAnalysis <- function(
	x = "ASEset",
	...
	){
	object <- new("GlobalAnalysis", data = list())

		if(!class(x)=="ASEset"){
			stop("x must be of class ASEset")
		}

	
	#Return object
	object
}

#' Initialize riskVariant
#' 
#' Functions to construct riskVariant objects
#' 
#' produces a class container for reference bias calculations
#' 
#' @name initialize-riskVariant
#' @rdname initialize-riskVariant
#' @aliases initialize-riskVariant 
#' @param x GRanges object for the SNPs
#' @param genotype matrix
#' @param colData A \code{DataFrame} object containing sample specific data
#' @param ... internal arguments
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords bias mapbias refBias
#' @examples
#'
#' data(ASEset)
#' ge <- inferGenotypes(ASEset)
#' rv <- riskVariantFromGRanges(x=GRvariants, genotype=ge)
#'
#' 
NULL

#' @rdname initialize-riskVariant
#' @export 
#setMethod("riskVariant","riskVariant", function(
riskVariantFromGRanges <- function(
	x,
	genotype,
	colData = NULL,
	...){

	if(is.null(colData)){
		sset <- SummarizedExperiment(
					assays = SimpleList(
						genotype=genotype
						), 
					rowRanges = x
				) 
	}else{
		sset <- SummarizedExperiment(
					assays = SimpleList(
						genotype=genotype
						), 
					rowRanges = x,
					colData = colData
				) 
	}

	rownames(sset) <- names(x)

	#valid
	#validObject(.Object)

	#Return object
	new("riskVariant", sset,
		meta = list()
	)
}


