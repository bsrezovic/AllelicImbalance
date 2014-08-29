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
#' Four different alleleCount options are available. The simples one is the
#' nonStranded option, and is experiments where the strand information is not
#' known. In this option both the plus, minus and unknown strand will be
#' counted and present. The unknown strand is when the aligner could not find
#' any strand associated with the read. Then there are an option too add plus
#' and minus stranded data. When using this, it is essential to make sure that
#' the RNA-seq experiment under analysis has in fact been created so that
#' correct strand information was obtained.
#' 
#' @name ASEset-class
#' @rdname ASEset-class
#' @aliases ASEset-class ASEset alleleCounts mapBias fraction arank
#' alleleCounts,ASEset-method mapBias,ASEset-method fraction,ASEset-method
#' arank,ASEset-method
#' @docType class
#' @param x ASEset object
#' @param strand which strand of 'nonStranded', '+', '-' or '*'
#' @param verbose makes function more talkative
#' @param ret return names or counts
#' @param ... additional arguments
#' @return An object of class ASEset containing location information and allele
#' counts for a number of SNPs measured in a number of samples on various
#' strand, as well as mapBias information. All data is stored in a manner
#' similar to the \code{\link[GenomicRanges]{SummarizedExperiment}} class.
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
#' \t#make example countList
#' \tset.seed(42)
#' \tcountListPlus <- list()
#' \tsnps <- c('snp1','snp2','snp3','snp4','snp5')
#' \tfor(snp in snps){
#' \t  count<-matrix(rep(0,16),ncol=4,dimnames=list(
#' \t\tc('sample1','sample2','sample3','sample4'),
#' \t\tc('A','T','G','C')))
#' \t  #insert random counts in two of the alleles 
#' \t  for(allele in sample(c('A','T','G','C'),2)){
#' \t\tcount[,allele]<-as.integer(rnorm(4,mean=50,sd=10))
#' \t  }
#' \t  countListPlus[[snp]] <- count
#' \t}
#' 
#' \t#make example rowData
#' \trowData <- GRanges(
#' \t  seqnames = Rle(c('chr1', 'chr2', 'chr1', 'chr3', 'chr1')),
#' \t  ranges = IRanges(1:5, width = 1, names = head(letters,5)),
#' \t  snp = paste('snp',1:5,sep='')
#' \t)
#' \t#make example colData
#' \tcolData <- DataFrame(Treatment=c('ChIP', 'Input','Input','ChIP'), 
#' \t \trow.names=c('ind1','ind2','ind3','ind4'))
#' 
#' \t#make ASEset 
#' \ta <- ASEsetFromCountList(rowData, countListPlus=countListPlus, 
#' \tcolData=colData)
#' 
#' 
#'
#' @exportClass ASEset
#' @exportMethod alleleCounts mapBias fraction arank
setClass("ASEset", contains = "SummarizedExperiment", representation(variants = "vector"))

#' @rdname ASEset-class
setGeneric("alleleCounts", function(x, strand = "nonStranded") {
    standardGeneric("alleleCounts")
})

setMethod("alleleCounts", signature(x = "ASEset"), function(x, strand = "nonStranded") {
    if (!sum(strand %in% c("+", "-", "*", "nonStranded")) > 0) {
        stop("strand parameter has to be either '+', '-', '*' or 'nonStranded' ")
    }
    
    if (strand == "+") {
        el <- "countsPlus"
    } else if (strand == "-") {
        el <- "countsMinus"
    } else if (strand == "*") {
        el <- "countsUnknown"
    } else if (strand == "nonStranded") {
        el <- "countsNonStranded"
    } else {
        stop("unknown strand option")
    }
    
    # check if strand option is present as assay
    if (!(el %in% names(assays(x)))) {
        stop("strand is not present as assay in ASEset object")
    }
    
    # assume alleleCount information is stored as element 1
    alleleCountList <- list()
    
    for (i in 1:nrow(assays(x)[[el]])) {
        mat <- assays(x)[[el]][i, , ]
        if (class(mat) == "integer") {
            mat <- t(as.matrix(mat))
            # rownames(mat) <- colnames(x)
        }
        if (class(mat) == "numeric") {
            mat <- t(mat)
            colnames(mat) <- x@variants
        } else {
            colnames(mat) <- x@variants
        }
        rownames(mat) <- colnames(x)
        alleleCountList[[i]] <- mat
    }
    # add snp id
    names(alleleCountList) <- rownames(x)
    
    # return object
    alleleCountList
    
})

#' @rdname ASEset-class
setGeneric("mapBias", function(x) {
    standardGeneric("mapBias")
})

setMethod("mapBias", signature(x = "ASEset"), function(x) {
    # assume alleleCount information is stored as element 1
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
    
    
    # return object
    mapBiasList
    
})

#' @rdname ASEset-class
setGeneric("fraction", function(x, strand = "nonStranded", verbose = FALSE) {
    standardGeneric("fraction")
})

setMethod("fraction", signature(x = "ASEset"), function(x, strand = "nonStranded", 
    verbose = FALSE) {
    
    if (!sum(strand %in% c("+", "-", "*", "nonStranded")) > 0) {
        stop("strand parameter has to be either '+', '-', '*' or 'nonStranded' ")
    }
    
    if (strand == "+") {
        el <- "countsPlus"
    } else if (strand == "-") {
        el <- "countsMinus"
    } else if (strand == "*") {
        el <- "countsUnknown"
    } else if (strand == "nonStranded") {
        el <- "countsNonStranded"
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
            cat(paste("Snp", "was possible tri-allelic, but only two most frequent alleles were plotted. Counts:"), 
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
setGeneric("arank", function(x, ret = "names", strand = "nonStranded", ...) {
    standardGeneric("arank")
})

setMethod("arank", signature(x = "ASEset"), function(x, ret = "names", strand = "nonStranded", 
    ...) {
    acounts <- alleleCounts(x, strand = strand)
    
    if (ret == "names") {
        arank <- lapply(acounts, function(x) {
            names(sort(apply(x, 2, sum), decreasing = TRUE))
        })
    }
    if (ret == "counts") {
        arank <- lapply(acounts, function(x) {
            as.numeric(sort(apply(x, 2, sum), decreasing = TRUE))
        })
    }
    arank
}) 
