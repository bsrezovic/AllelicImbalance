#'@include chisq.test-methods.R
NULL

#' binomial test
#' 
#' Performs a binomial test on an ASEset object.
#' 
#' the test can only be applied to one strand at the time.
#' 
#' @name binom.test
#' @rdname binom.test
#' @aliases binom.test,ASEset-method
#' @docType methods
#' @param x \code{ASEset} object
#' @param n strand option
#' @return \code{binom.test} returns a matrix
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{chisq.test}} which is another test
#' that can be applied on an \link{ASEset} object.  }
#' @keywords binomial test
#' @examples
#' 
#' \t#load example data
#' \tdata(ASEset)
#' \t
#' \t#make a binomial test
#' \tbinom.test(ASEset,'nonStranded')
#' 
#' 
#'  @importFrom stats binom.test
#'  @exportMethod binom.test
NULL

#' @rdname binom.test
setGeneric("binom.test")

#' @rdname binom.test
setMethod("binom.test", signature(x = "ASEset", n = "ANY"), function(x, n = "nonStranded") {
    
    strand <- n
    
    # checks
    if (!sum(strand %in% c("+", "-", "*", "nonStranded")) > 0) {
        stop("strand parameter n has to be either '+', '-', '*' or 'nonStranded' ")
    }
    
    
    biasWarning <- vector()
    pLst <- list()
    
    # use this strand to get dimension information
    if (strand == "both") {
        tmpStrand <- "+"
    } else {
        tmpStrand <- strand
    }
    
    for (i in 1:length(alleleCounts(x, strand = tmpStrand))) {
        bias <- mapBias(x)[[i, drop = FALSE]]
        
        df <- alleleCounts(x, tmpStrand)[[i]]
        returnVec <- rep(NA, nrow(df))
        for (j in 1:nrow(df)) {
            so <- sort(df[j, ], decreasing = TRUE)
            if (so[2] != 0) {
                # place bias in same order as so
                bi <- bias[, names(so)[1:2]]
                
                # check that bi[1] and bi[2] adds up to 100% or force 0.5-0.5 and give warning
                if (bi[1] + bi[2] != 1) {
                  bi[1] <- bi[2] <- 0.5
                  warning("Found disrepancy between mapping bias information and mapBiasExpMean and the allele count. Coerced to expectation of 0.5 to 0.5 ratio")
                }
                expCounts <- c(sum(so) * bi[1], sum(so) * bi[2])
                returnVec[j] <- binom.test(as.numeric(c(so[1], so[2])), p = as.numeric(bi[1]/(bi[1] + 
                  bi[2])))[[3]]
            } else {
                returnVec[j] <- NA
            }  # if it's mono-allelic
        }  #end for-loop
        pLst[[names(alleleCounts(x, strand = tmpStrand)[i])]] <- returnVec
    }
    if (length(biasWarning) > 0) {
        warning(paste(length(biasWarning), "SNPs had disrepancy between their two highest read count alleles and the\n\n\t\t\t\t\t\t\t\t\ttwo alleles indicated in mapBiasExpMean. They were coerced to a standard and expected alignment \n\t\t\t\t\t\t\t\t\tbias of 0.5 to 0.5, but could be checked further as this indicates non-standard allele distributions:\n", 
            paste(biasWarning[1:(min(c(length(biasWarning), 20)))], collapse = ", ")))
    }
    return(as.matrix(as.data.frame(pLst)))
}) 
