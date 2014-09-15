#'@include auxillary-functions.R
NULL

#' chi-square test
#' 
#' Performs a chisq.test on an ASEset object.
#' 
#' The test is performed on one strand in an ASEset object.
#' 
#' @name chisq.test
#' @rdname chisq.test
#' @aliases chisq.test,ASEset-method
#' @docType methods
#' @param x \code{ASEset} object
#' @param y strand option
#' @param p NOT USED
#' @param correct NOT USED
#' @param rescale.p NOT USED
#' @param simulate.p.value NOT USED
#' @param B NOT USED
#' @return \code{chisq.test} returns a matrix with the chisq.test P-value for
#' each SNP and sample
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{binom.test}} which is another test
#' that can be applied on an \link{ASEset} object.  }
#' @keywords chi-square test
#' @examples
#' 
#' #load example data
#' data(ASEset)
#' 
#' #make a chi-square test on default non-stranded strand 
#' chisq.test(ASEset)
#' 
#'  @importFrom stats chisq.test
#'  @exportMethod chisq.test
NULL

#' @rdname chisq.test
setGeneric("chisq.test")

#' @rdname chisq.test
setMethod("chisq.test", signature(x = "ASEset", y = "ANY"), function(x, y = "*") {
    strand = y
    
    if (!sum(strand %in% c("+", "-", "*")) > 0) {
        stop(paste("strand parameter y has to be either '+', '-' or '*'"))
    }
    
    biasWarning <- vector()
    pLst <- list()
    
    if (strand == "both") {
        tmpStrand <- "+"
    } else {
        tmpStrand <- strand
    }
    
    for (i in 1:length(alleleCounts(x, strand = tmpStrand))) {
        
        bias <- mapBias(x)[[i, drop = FALSE]]
        
        df <- alleleCounts(x, strand = tmpStrand)[[i]]
        
        returnVec <- rep(NA, nrow(df))
        for (j in 1:nrow(df)) {
            so <- sort(df[j, ], decreasing = TRUE)
            if (so[2] != 0) {
                # place bias in same order as so
                bi <- bias[j, names(so)[1:2]]
                
                # check that bi[1] and bi[2] adds up to 100% or force 0.5-0.5 and give warning
                if (bi[1] + bi[2] != 1) {
                  bi[1] <- bi[2] <- 0.5
                  warning("Found disrepancy between mapping bias information and mapBiasExpMean and the allele count. Coerced to expectation of 0.5 to 0.5 ratio")
                }
                expCounts <- c(sum(so) * bi[1], sum(so) * bi[2])
                if (!(expCounts[1] < 5 | expCounts[2] < 5)) {
                  # the condition is that none of the Expected Counts should be below 5
                  returnVec[j] <- chisq.test(c(so[1], so[2]), p = c(bi[1], bi[2]), 
                    correct = FALSE)[[3]]
                } else {
                  returnVec[j] <- NA
                }  #if it's below 5 expected counts in either 
            } else {
                returnVec[j] <- NA
            }  # if it's mono-allelic
        }  #end for-loop
        pLst[[rownames(x)[i]]] <- returnVec
    }
    if (length(biasWarning) > 0) {
        warning(paste(length(biasWarning), "SNPs had disrepancy between their two highest read count alleles and the\n\n\t\t\t\t\t\t\t\t\ttwo alleles indicated in mapBiasExpMean. They were coerced to a standard and expected alignment \n\t\t\t\t\t\t\t\t\tbias of 0.5 to 0.5, but could be checked further as this indicates non-standard allele distributions:\n", 
            paste(biasWarning[1:(min(c(length(biasWarning), 20)))], collapse = ", ")))
    }
    return(as.matrix(as.data.frame(pLst)))
}) 
