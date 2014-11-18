#'@include AllelicImbalance-package.R
NULL

# setMethods('initialize') is not required at the moment

#' Initialize ASEset
#' 
#' Functions to construct ASEset objects
#' 
#' The resulting \code{ASEset} object is based on the
#' \code{SummarizedExperiment}, and will therefore inherit the same accessors
#' and ranges operations.
#' 
#' If both countListPlus and countListMinus are given they will be used to 
#' calculate countListUnknown, which is the sum of the plus and minus strands.
#' 
#' countListPlus, countListMinus and countListUnknown are
#' i.e. the outputs from the \code{\link{getAlleleCounts}} function.
#' 
#' @name initialize-ASEset
#' @rdname initialize-ASEset
#' @aliases initialize-ASEset ASEsetFromCountList
#' @param rowData A \code{GenomicRanges object} that contains the variants of
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
#' @param colData A \code{DataFrame} object containing sample specific data
#' @param mapBiasExpMean A 3D \code{array} where the SNPs are in the 1st
#' dimension, samples in the 2nd dimension and variants in the 3rd dimension.
#' @param verbose Makes function more talkative
#' @param ... arguments passed on to SummarizedExperiment constructor
#' @return \code{ASEsetFromCountList} returns an \code{ASEset} object.
#' @note \code{ASEsetFromCountList} requires the same input data as an
#' SummarizedExperiment, but with minimum one assay for the allele counts.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The
#' \code{\link[GenomicRanges]{SummarizedExperiment}} for ranges operations.  }
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
#' #make example rowData
#' rowData <- GRanges(
#' seqnames = Rle(c('chr1', 'chr2', 'chr1', 'chr3', 'chr1')),
#'          ranges = IRanges(1:5, width = 1, names = head(letters,5)),
#'          snp = paste('snp',1:5,sep='')
#'          )
#' #make example colData
#' colData <- DataFrame(Treatment=c('ChIP', 'Input','Input','ChIP'), 
#'  row.names=c('ind1','ind2','ind3','ind4'))
#' 
#' #make ASEset 
#' a <- ASEsetFromCountList(rowData, countListPlus=countListPlus, 
#' countListMinus=countListMinus, colData=colData)
#' 
NULL


#' @rdname initialize-ASEset
#' @export 
ASEsetFromCountList <- function(rowData, countListUnknown = NULL, countListPlus = NULL, 
    countListMinus = NULL, colData = NULL, mapBiasExpMean = NULL, 
    verbose = FALSE, ...) {
    
    if (verbose) {
        cat("rowData\n")
        cat(class(rowData))
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
    m <- matrix(unlist(l), ncol = length(l))
    if (sum(!(m == m[1])) > 0) {
        stop("one or more list contains more or less columns than the others")
    }
    
    # check that the number of rows in each dataframe are the same
    l <- unlist(lapply(countLists, function(x) {
        lapply(get(x), nrow)
    }))
    m <- matrix(unlist(l), ncol = length(l))
    if (sum(!(m == m[1])) > 0) {
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
        m <- matrix(unlist(l), ncol = length(l))
        if (sum(!(m == m[, 1])) > 0) {
            stop(paste("the names or order of names are not the same in all data frames in list", 
                countListName))
        }
        cMatrix[, i] <- m[, 1]
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
        m <- matrix(unlist(l), ncol = length(l))
        if (sum(!(m == m[, 1])) > 0) {
            stop(paste("the names or order of names are not the same in all data frames in list", 
                countListName))
        }
        cMatrix[, i] <- m[, 1]
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
        
        for (i in 1:(dim(m)[2])) {
            # check that no snp and sample exceeds sum 1 frequency
            if (!sum(!apply(m[, i, ], 1, sum) == 1) == 0) {
                stop(paste("for each snp and sample the sum of allele frequencies must sum to one,\n\t\t\t\tfor sample element", 
                  i, "and snp nr", paste(which(!apply(m[, i, ], 1, sum) == 1), collapse = " "), 
                  "this was not fullfilled"))
            }
            # check for tri-allelic cases, which we dont allow as mapping biases.
            if (!sum(!(apply((m[, i, ] > 0), 1, sum) == 2)) == 0) {
                stop(paste("tri-allelic SNPs have not been implemented yet. Please write an email if this is of interest.\n \n\t\t\t\tTri-allelic case was found for sample nr", 
                  i, "and snp nr", paste(which(!apply((m[, i, ] > 0), 1, sum) == 
                    2), collapse = " ")))
            }
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
        assays[["countsPlus"]] <- ar1
    }
    # minus
    if (!is.null(countListMinus)) {
		#empty array that handles only four nucleotides 
        ar2 <- array(NA, c(snps, ind, 4))  
        for (i in 1:snps) {
            ar2[i, , ] <- countListMinus[[i]]
        }
        assays[["countsMinus"]] <- ar2
        
    }
    # unknown
    if (!is.null(countListUnknown)) {
		#empty array that handles only four nucleotides 
        ar3 <- array(NA, c(snps, ind, 4))  
        for (i in 1:snps) {
            ar3[i, , ] <- countListUnknown[[i]]
        }
        assays[["countsUnknown"]] <- ar3
        
    }else if((!is.null(countListMinus)) & (!is.null(countListPlus))){
		#Calculate the non-stranded representative from the stranded data
        assays[["countsUnknown"]] <- ar1+ar2
	}
    
    # assign mapBiasExpMean
    if (is.null(mapBiasExpMean)) {
        assays[["mapBias"]] <- getDefaultMapBiasExpMean3D(countList)
    } else {
        assays[["mapBias"]] <- mapBiasExpMean
    }
    
    if (is.null(colData)) {
        colData <- DataFrame(row.names = unlist(unique(lapply(countList, rownames))))
    }
    
    sset <- SummarizedExperiment(assays = assays, rowData = rowData, colData = colData, 
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
ASEsetFromArrays <- function(rowData, countsUnknown = NULL, countsPlus = NULL, 
    countsMinus = NULL, colData = NULL, mapBiasExpMean = NULL, 
    verbose = FALSE, ...) {
    
    if (verbose) {
        cat("rowData\n")
        cat(class(rowData))
        cat("countsPlus\n")
        cat(class(countsPlus))
        cat("countsMinus\n")
        cat(class(countsMinus))
        cat("countsUnknown\n")
        cat(class(countsUnknown))
        cat("colData\n")
        cat(class(colData))
        cat("mapBiasExpMean\n")
        cat(class(mapBiasExpMean))
    }
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
        
        for (i in 1:(dim(m)[2])) {
            # check that no snp and sample exceeds sum 1 frequency
            if (!sum(!apply(m[, i, ], 1, sum) == 1) == 0) {
                stop(paste("for each snp and sample the sum of allele frequencies must sum to one,\n\t\t\t\tfor sample element", 
                  i, "and snp nr", paste(which(!apply(m[, i, ], 1, sum) == 1), collapse = " "), 
                  "this was not fullfilled"))
            }
            # check for tri-allelic cases, which we dont allow as mapping biases.
            if (!sum(!(apply((m[, i, ] > 0), 1, sum) == 2)) == 0) {
                stop(paste("tri-allelic SNPs have not been implemented yet. Please write an email if this is of interest.\n \n\t\t\t\tTri-allelic case was found for sample nr", 
                  i, "and snp nr", paste(which(!apply((m[, i, ] > 0), 1, sum) == 
                    2), collapse = " ")))
            }
        }
    }
    
    # choose a common countList by picking the first one, for dimension info
    countList <- get(countLists[1])
    ind <- length(unlist(unique(lapply(countList, rownames))))
    snps <- length(countList)
    
    # SimpleList init
    assays <- SimpleList()
    
    # plus
    if (!is.null(countsPlus)) {
		#empty array that handles only four nucleotides 
        #ar1 <- array(NA, c(snps, ind, 4))  
        #for (i in 1:snps) {
        #    ar1[i, , ] <- countsPlus[[i]]
        #}
        assays[["countsPlus"]] <- countsPlus
    }
    # minus
    if (!is.null(countsMinus)) {
        assays[["countsMinus"]] <- countsMinus
        
    }
    # unknown
    if (!is.null(countsUnknown)) {
        assays[["countsUnknown"]] <- countsUnknown
        
    }else if((!is.null(countsMinus)) & (!is.null(countsPlus))){
		#Calculate the non-stranded representative from the stranded data
        assays[["countsUnknown"]] <- countsPlus+countsMinus
	}
    
    # assign mapBiasExpMean
    if (is.null(mapBiasExpMean)) {
        assays[["mapBias"]] <- getDefaultMapBiasExpMean3D(countList)
    } else {
        assays[["mapBias"]] <- mapBiasExpMean
    }
    
    if (is.null(colData)) {
        colData <- DataFrame(row.names = unlist(unique(lapply(countList, rownames))))
    }
    
    sset <- SummarizedExperiment(assays = assays, rowData = rowData, colData = colData, 
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

	sset <- SummarizedExperiment(assays = SimpleList(referenceFrequency=assay), rowData = rowData(x), colData = colData(x)) 
	rownames(sset) <- rownames(x)

	#valid
	#validObject(.Object)

	#Return object
	new("ReferenceBias", sset, strands = c("*","+","-"))
}


##' @rdname initialize-ReferenceBias
#refBias <- function(x){
#
#		if(!class(x)=="ASEset"){
#			stop("x must be of class ASEset")
#		}
#
#        # create object
#        
#}
#

