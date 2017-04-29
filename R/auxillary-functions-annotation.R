#'@include ASEset-class.R
NULL

#' AnnotationDb wrappers
#' 
#' These functions acts as wrappers to retrieve information from annotation
#' database objects (\code{annotationDb objects}) or (\code{transcriptDb
#' objects})
#' 
#' These functions retrieve regional annotation from OrgDb or TxDb objects,
#' when given GRanges objects.
#' 
#' @name annotation-wrappers
#' @rdname annotation-wrappers
#' @aliases getGenesFromAnnotation getGenesVector getExonsFromAnnotation
#' getTranscriptsFromAnnotation getCDSFromAnnotation getExonsVector
#' getTranscriptsVector getCDSVector getAnnotationDataFrame
#' @param OrgDb An \code{OrgDb} object
#' @param GR A \code{GenomicRanges} object with sample area
#' @param leftFlank An \code{integer} specifying number of additional
#' nucleotides around the SNPs for the leftFlank
#' @param rightFlank An \code{integer} specifying number of additional
#' nucleotides around the SNPs for the rightFlank
#' @param getUCSC A \code{logical} indicating if UCSC transcript IDs should
#' also be retrieved
#' @param TxDb A \code{transcriptDb} object
#' @param strand Two options,'+' or '-'
#' @param annotationType select one or more from 'gene', 'exon', 'transcript',
#' 'cds'.
#' @param verbose A \code{logical} making the functions more talkative
#' @return %OrgDb The \code{getGenesFromAnnotation} function will return a
#' \code{GRanges object} with ranges over the genes in the region.
#' 
#' The \code{getGenesVector} function will return a character vector where each
#' element are gene symbols separated by comma
#' 
#' %transcriptDb The \code{getExonsFromAnnotation} function will return a
#' \code{GRanges object} with ranges over the exons in the region.
#' 
#' The \code{getTranscriptsFromAnnotation} function will return a \code{GRanges
#' object} with ranges over the transcripts in the region.
#' 
#' The \code{getCDSFromAnnotation} function will return a \code{GRanges object}
#' with ranges over the CDSFs in the region.
#' 
#' The \code{getExonsVector} function will return a character vector where each
#' element are exons separated by comma
#' 
#' The \code{getTranscriptsVector} function will return a character vector
#' where each element are transcripts separated by comma
#' 
#' The \code{getCDSVector} function will return a character vector where each
#' element are CDSs separated by comma
#' 
#' The \code{getAnnotationDataFrame} function will return a data.frame with
#' annotations. This function is used internally by i.e. the barplot-function
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords genes exons transcripts CDS annotation
#' @examples
#' 
#' 
#'   data(ASEset)
#'   require(org.Hs.eg.db)
#'   require(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'   OrgDb <- org.Hs.eg.db
#'   TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' 
#'   #use for example BcfFiles as the source for SNPs of interest
#'   GR <- rowRanges(ASEset)
#'   #get annotation
#'   g <- getGenesFromAnnotation(OrgDb,GR)
#'   e <- getExonsFromAnnotation(TxDb,GR)
#'   t <- getTranscriptsFromAnnotation(TxDb,GR)
#'   c <- getCDSFromAnnotation(TxDb,GR)
#'   
#' @export getGenesFromAnnotation
#' @export getExonsFromAnnotation
#' @export getTranscriptsFromAnnotation
#' @export getCDSFromAnnotation
#' @export getGenesVector
#' @export getExonsVector
#' @export getTranscriptsVector
#' @export getCDSVector
#' @export getAnnotationDataFrame
NULL

#' @rdname annotation-wrappers
getGenesFromAnnotation <- function(OrgDb, GR, leftFlank = 0, rightFlank = 0, getUCSC = FALSE, 
    verbose = FALSE) {
    # checks
    if (class(OrgDb) != "OrgDb") 
        stop(paste("OrgDb must of class OrgDb, not", class(OrgDb)))
    
    if (class(GR) != "GRanges") 
        stop(paste("GR must of class GRanges, not", class(GR)))
    
    if (!class(leftFlank) %in% c("numeric")) 
        stop(paste("leftFlank must be of class numeric, not:", class(leftFlank)))
    if (length(leftFlank) != 1) 
        stop(paste("leftFlank must be of length 1, not:", length(leftFlank)))
    if (leftFlank < 0) 
        stop(paste("leftFlank must be equal to or larger than 0"))
    
    if (!class(rightFlank) %in% c("numeric")) 
        stop(paste("rightFlank must be of class numeric, not:", class(rightFlank)))
    if (length(rightFlank) != 1) 
        stop(paste("rightFlank must be of length 1, not:", length(rightFlank)))
    if (rightFlank < 0) 
        stop(paste("rightFlank must be equal to or larger than 0"))
    
    if (!class(getUCSC) %in% c("logical")) 
        stop(paste("getUCSC must be of class logical, not:", class(getUCSC)))
    if (length(getUCSC) != 1) 
        stop(paste("getUCSC must be of length 1, not:", length(getUCSC)))
    
    if (!"UCSCKG" %in% columns(OrgDb)) {
        if (verbose) 
            message("Unable to retrieve UCSCKG column from OrgDb. Omitting")
        getUCSC <- FALSE
    }
    
    if (!class(verbose) %in% c("logical")) 
        stop(paste("verbose must be of class logical, not:", class(verbose)))
    if (length(verbose) != 1) 
        stop(paste("verbose must be of length 1, not:", length(verbose)))
    
    # remove chr in seqnames
    seqLevels <- sub("^chr", "", seqlevels(GR))
    seqlevels(GR) <- seqLevels
    
    # pre-filtering to local region +/- 1MB (for speed purposes)
    startFilter <- max(c(1, start(range(GR)) - 10^6))
    endFilter <- end(range(GR)) + 10^6
    colsFilter <- c("CHR", "CHRLOC", "CHRLOCEND", "SYMBOL")
    sFilter <- suppressWarnings(select(OrgDb, keys = seqLevels, columns = colsFilter, 
        keytype = "CHR"))
    symbolsToGet <- sFilter[abs(sFilter[, "CHRLOC"]) > startFilter & abs(sFilter[, 
        "CHRLOCEND"]) < endFilter & !is.na(sFilter[, "CHRLOCEND"]) & !is.na(sFilter[, 
        "CHRLOC"]), "SYMBOL"]
    
    
    # then create an annGR for genes in range
    if (getUCSC) {
        cols <- c("SYMBOL", "CHR", "CHRLOC", "CHRLOCEND", "ENSEMBL", "UCSCKG")
    } else {
        cols <- c("SYMBOL", "CHR", "CHRLOC", "CHRLOCEND", "ENSEMBL")
    }
    s <- suppressWarnings(select(OrgDb, keys = symbolsToGet, columns = cols, keytype = "SYMBOL"))
    
    
    # remove Symbols with NAs
    TFminusStrand <- s[["CHRLOC"]] < 0
    TFplusStrand <- s[["CHRLOC"]] > 0
    sNoNas <- s[c(which(TFminusStrand), which(TFplusStrand)), ]
    
    # make Strand vector
    TFminusStrand2 <- sNoNas[["CHRLOC"]] < 0
    strand <- rep("+", length = (dim(sNoNas)[1]))
    
    strand[TFminusStrand2] <- "-"
    
    # make start and end vector
    sNonNegative <- sNoNas
    sNonNegative[TFminusStrand2, c("CHRLOC", "CHRLOCEND")] <- -sNonNegative[TFminusStrand2, 
        c("CHRLOC", "CHRLOCEND")]
    start <- sNonNegative[["CHRLOC"]]
    end <- sNonNegative[["CHRLOCEND"]]
    
    # make seqnames
    seqnames <- sNonNegative[["CHR"]]
    
    # make the annGR containing all genes
    if (getUCSC) {
        annGR <- GRanges(seqnames = Rle(seqnames), ranges = IRanges(start, end), 
            strand = Rle(strand), Symbol = sNonNegative[["SYMBOL"]], Ensembl = sNonNegative[["ENSEMBL"]], 
            UCSCKG = sNonNegative[["UCSCKG"]])
    } else {
        annGR <- GRanges(seqnames = Rle(seqnames), ranges = IRanges(start, end), 
            strand = Rle(strand), Symbol = sNonNegative[["SYMBOL"]], Ensembl = sNonNegative[["ENSEMBL"]])
        
    }
    
    # check that all levels in GR exist in annGR, if not exclude these levels
    if (sum(!levels(seqnames(GR)) %in% levels(seqnames(annGR))) > 0) {
        TFkeepLevels <- levels(seqnames(GR)) %in% levels(seqnames(annGR))
        seqlevels(GR) <- seqlevels(GR)[TFkeepLevels]
        
    }
    
    # check that all levels in annGR exist in GR, if not exclude these levels
    if (sum(!levels(seqnames(annGR)) %in% levels(seqnames(GR))) > 0) {
        TFkeepLevels <- levels(seqnames(annGR)) %in% levels(seqnames(GR))
        seqlevels(annGR, pruning.mode="coarse") <- seqlevels(annGR)[TFkeepLevels]
    }
    
    # the seqlevels comes in different orders. This will give the correct order.
    seqlevels(GR) <- seqlevels(annGR)
    
    # find overlaps between annGR and Snps incl. flank region
    GenesInRegion <- subsetByOverlaps(annGR, GR)  #put in flankSize here when you have time ;)
    seqlengths(GenesInRegion) <- seqlengths(GR)
    GenesInRegion
}

#' @rdname annotation-wrappers
getGenesVector <- function(OrgDb, GR, leftFlank = 0, rightFlank = 0, verbose = FALSE) {
    if (verbose) {
        cat("start gene extraction\n")
    }
    GenesInRegion <- getGenesFromAnnotation(OrgDb, GR, leftFlank = leftFlank, rightFlank = rightFlank, 
        verbose = FALSE)
    
    seqlevels(GR) <- seqlevels(GenesInRegion)
    
    # remove duplicate symbol names if same symbol merge regions
    symbolList <- unique(mcols(GenesInRegion)[["Symbol"]])
    newGenesInRegion <- GRanges()
    
    # check if GenesInRegion is zero
    if (!length(GenesInRegion) == 0) {
        for (i in 1:length(symbolList)) {
            symbol <- symbolList[i]
            TF <- mcols(GenesInRegion)[["Symbol"]] == symbol
            
            G <- GRanges(seqnames = unique(seqnames(GenesInRegion[TF])), ranges = IRanges(min(start(GenesInRegion[TF])), 
                max(end(GenesInRegion[TF]))), strand = unique(strand(GenesInRegion[TF])))
            mcols(G) <- unique(mcols(GenesInRegion[TF])[, "Symbol", drop = FALSE])
            
            newGenesInRegion <- c(newGenesInRegion, G)
        }
    }
    
    # half-vectorized solution
    h <- findOverlaps(newGenesInRegion, GR)
    symbolVec <- vector()
    for (i in 1:length(GR)) {
        symbolVec[i] <- paste(mcols(newGenesInRegion[queryHits(h)[subjectHits(h) == 
            i]])[["Symbol"]], collapse = ",")
    }
    # set NAs where appropriate
    symbolVec[symbolVec == ""] <- NA
    # return list with symbols
    symbolVec
    
    # return list with symbols
    symbolVec
}

#' @rdname annotation-wrappers
getExonsFromAnnotation <- function(TxDb, GR, leftFlank = 0, rightFlank = 0, verbose = FALSE) {
    
    # checks
    if (class(TxDb) != "TxDb") 
        stop(paste("GR must of class TxDb, not", class(TxDb)))
    
    if (class(GR) != "GRanges") 
        stop(paste("GR must of class GRanges, not", class(GR)))
    
    if (!class(leftFlank) %in% c("numeric")) 
        stop(paste("leftFlank must be of class numeric, not:", class(leftFlank)))
    if (length(leftFlank) != 1) 
        stop(paste("leftFlank must be of length 1, not:", length(leftFlank)))
    if (leftFlank < 0) 
        stop(paste("leftFlank must be equal to or larger than 0"))
    
    if (!class(rightFlank) %in% c("numeric")) 
        stop(paste("rightFlank must be of class numeric, not:", class(rightFlank)))
    if (length(rightFlank) != 1) 
        stop(paste("rightFlank must be of length 1, not:", length(rightFlank)))
    if (rightFlank < 0) 
        stop(paste("rightFlank must be equal to or larger than 0"))
    
    if (!class(verbose) %in% c("logical")) 
        stop(paste("verbose must be of class logical, not:", class(verbose)))
    if (length(verbose) != 1) 
        stop(paste("verbose must be of length 1, not:", length(verbose)))
    
    # remove chr in seqnames for GR
    seqLevels <- sub("^chr", "", seqlevels(GR))
    seqlevels(GR) <- seqLevels
    
    seqlevels(TxDb, pruning.mode="coarse") <- paste("chr", names(seqlengths(GR)), sep = "")
    
    # Get all exons from the active chromosomes By creating a GRanges from TxDb
    annGR <- exons(TxDb, columns = c("exon_id", "tx_name"))
    
    # remove chr in seqnames for annGR
    seqLevels <- sub("^chr", "", seqlevels(annGR))
    seqlevels(annGR) <- seqLevels
    
    # check that all levels in GR exist in annGR, if not exclude these levels
    if (sum(!levels(seqnames(GR)) %in% levels(seqnames(annGR))) > 0) {
        TFkeepLevels <- levels(seqnames(GR)) %in% levels(seqnames(annGR))
        seqlevels(GR) <- seqlevels(GR)[TFkeepLevels]
        
    }
    
    # check that all levels in annGR exist in GR, if not exclude these levels
    if (sum(!levels(seqnames(annGR)) %in% levels(seqnames(GR))) > 0) {
        TFkeepLevels <- levels(seqnames(annGR)) %in% levels(seqnames(GR))
        seqlevels(annGR, pruning.mode="coarse") <- seqlevels(annGR)[TFkeepLevels]
    }
    
    # the seqlevels comes in different orders. This will give the correct order.
    seqlevels(GR) <- seqlevels(annGR)
    
    # add flanking regions
    lf <- flank(GR, leftFlank, start = TRUE)
    rf <- flank(GR, rightFlank, start = FALSE)
    GR <- c(lf, GR, rf)
    start <- max(c(1, min(start(GR))))
    end <- max(end(GR))
    
    # speed-increasing coarse subset to +/- 1MB before extracting (because smaller
    # annGR is faster and we have to access several times)
    annGR <- annGR[start(ranges(annGR)) > max(c(1, start - 10^6)) & end(ranges(annGR)) < 
        (end + 10^6)]
    
    # extract names of all transcripts with exons in plotting window
    tx_names <- unique(unlist(as.list(mcols(annGR[start(ranges(annGR)) > start & 
        end(ranges(annGR)) < end])[["tx_name"]])))
    
    # retrieve all transcript annotation, including potential off-plot tails and
    # heads
    ExonsInRegion <- annGR[sapply(mcols(annGR)[["tx_name"]], function(x) {
        any(tx_names %in% x)
    })]
    
    seqlengths(ExonsInRegion) <- seqlengths(GR)
    ExonsInRegion
    
}

#' @rdname annotation-wrappers
getExonsVector <- function(TxDb, GR, leftFlank = 0, rightFlank = 0, verbose = FALSE) {
    if (verbose) {
        cat("start exon extraction\n")
    }
    
    ExonsInRegion <- getExonsFromAnnotation(TxDb, GR, leftFlank, rightFlank, verbose = verbose)
    
    seqlevels(GR) <- seqlevels(ExonsInRegion)
    
    # half-vectorized solution
    h <- findOverlaps(ExonsInRegion, GR)
    ExonVec <- vector()
    for (i in 1:length(GR)) {
        ExonVec[i] <- paste(mcols(ExonsInRegion[queryHits(h)[subjectHits(h) == i]])[["exon_id"]], 
            collapse = ",")
    }
    # set NAs where appropriate
    ExonVec[ExonVec == ""] <- NA
    # return list with symbols
    ExonVec
    
}

#' @rdname annotation-wrappers
getTranscriptsFromAnnotation <- function(TxDb, GR, leftFlank = 0, rightFlank = 0, 
    verbose = FALSE) {
    
    # checks
    if (class(TxDb) != "TxDb") 
        stop(paste("GR must of class TxDb, not", class(TxDb)))
    
    if (class(GR) != "GRanges") 
        stop(paste("GR must of class GRanges, not", class(GR)))
    
    if (!class(leftFlank) %in% c("numeric")) 
        stop(paste("leftFlank must be of class numeric, not:", class(leftFlank)))
    if (length(leftFlank) != 1) 
        stop(paste("leftFlank must be of length 1, not:", length(leftFlank)))
    if (leftFlank < 0) 
        stop(paste("leftFlank must be equal to or larger than 0"))
    
    if (!class(rightFlank) %in% c("numeric")) 
        stop(paste("rightFlank must be of class numeric, not:", class(rightFlank)))
    if (length(rightFlank) != 1) 
        stop(paste("rightFlank must be of length 1, not:", length(rightFlank)))
    if (rightFlank < 0) 
        stop(paste("rightFlank must be equal to or larger than 0"))
    
    if (!class(verbose) %in% c("logical")) 
        stop(paste("verbose must be of class logical, not:", class(verbose)))
    if (length(verbose) != 1) 
        stop(paste("verbose must be of length 1, not:", length(verbose)))
    
    
    
    # remove chr in seqnames for GR
    seqLevels <- sub("^chr", "", seqlevels(GR))
    seqlevels(GR) <- seqLevels
    
    seqlevels(TxDb, pruning.mode="coarse") <- paste("chr", names(seqlengths(GR)), sep = "")
    
    
    # Get all exons from the active chromosomes By creating a GRanges from TxDb
    annGR <- transcripts(TxDb)
    
    
    # remove chr in seqnames for annGR
    seqLevels <- sub("^chr", "", seqlevels(annGR))
    seqlevels(annGR) <- seqLevels
    
    # check that all levels in GR exist in annGR, if not exclude these levels
    if (sum(!levels(seqnames(GR)) %in% levels(seqnames(annGR))) > 0) {
        TFkeepLevels <- levels(seqnames(GR)) %in% levels(seqnames(annGR))
        seqlevels(GR) <- seqlevels(GR)[TFkeepLevels]
        
    }
    
    # check that all levels in annGR exist in GR, if not exclude these levels
    if (sum(!levels(seqnames(annGR)) %in% levels(seqnames(GR))) > 0) {
        TFkeepLevels <- levels(seqnames(annGR)) %in% levels(seqnames(GR))
        seqlevels(annGR, pruning.mode="coarse") <- seqlevels(annGR)[TFkeepLevels]
    }
    
    # the seqlevels comes in different orders. This will give the correct order.
    seqlevels(GR) <- seqlevels(annGR)
    
    # add flanking regions
    lf <- flank(GR, leftFlank, start = TRUE)
    rf <- flank(GR, rightFlank, start = FALSE)
    GR <- c(lf, GR, rf)
    # find overlaps between annGR and Snps incl. flank region
    TxInRegion <- subsetByOverlaps(annGR, GR)  #put in flankSize here when you have time ;)
    seqlengths(TxInRegion) <- seqlengths(GR)
    TxInRegion
}

#' @rdname annotation-wrappers
getTranscriptsVector <- function(TxDb, GR, leftFlank = 0, rightFlank = 0, verbose = FALSE) {
    if (verbose) {
        cat("start transcript extraction\n")
    }
    
    TxInRegion <- getTranscriptsFromAnnotation(TxDb, GR, leftFlank, rightFlank)
    
    seqlevels(GR) <- seqlevels(TxInRegion)
    
    # half-vectorized solution
    h <- findOverlaps(TxInRegion, GR)
    TxVec <- vector()
    for (i in 1:length(GR)) {
        TxVec[i] <- paste(mcols(TxInRegion[queryHits(h)[subjectHits(h) == i]])[["tx_id"]], 
            collapse = ",")
    }
    # set NAs where appropriate
    TxVec[TxVec == ""] <- NA
    # return list with symbols
    TxVec
    
}

#' @rdname annotation-wrappers
getCDSFromAnnotation <- function(TxDb, GR, leftFlank = 0, rightFlank = 0, verbose = FALSE) {
    # CDS are the coding regions that do not only code for proteins, but other also
    # other types like RNA.
    
    # checks
    if (class(TxDb) != "TxDb") 
        stop(paste("GR must of class TxDb, not", class(TxDb)))
    
    if (class(GR) != "GRanges") 
        stop(paste("GR must of class GRanges, not", class(GR)))
    
    if (!class(leftFlank) %in% c("numeric")) 
        stop(paste("leftFlank must be of class numeric, not:", class(leftFlank)))
    if (length(leftFlank) != 1) 
        stop(paste("leftFlank must be of length 1, not:", length(leftFlank)))
    if (leftFlank < 0) 
        stop(paste("leftFlank must be equal to or larger than 0"))
    
    if (!class(rightFlank) %in% c("numeric")) 
        stop(paste("rightFlank must be of class numeric, not:", class(rightFlank)))
    if (length(rightFlank) != 1) 
        stop(paste("rightFlank must be of length 1, not:", length(rightFlank)))
    if (rightFlank < 0) 
        stop(paste("rightFlank must be equal to or larger than 0"))
    
    if (!class(verbose) %in% c("logical")) 
        stop(paste("verbose must be of class logical, not:", class(verbose)))
    if (length(verbose) != 1) 
        stop(paste("verbose must be of length 1, not:", length(verbose)))
    
    # remove chr in seqnames for GR
    seqLevels <- sub("^chr", "", seqlevels(GR))
    seqlevels(GR) <- seqLevels
    
    seqlevels(TxDb, pruning.mode="coarse") <- paste("chr", names(seqlengths(GR)), sep = "")
    
    
    # Get all exons from the active chromosomes By creating a GRanges from TxDb
    annGR <- cds(TxDb)
    
    
    # remove chr in seqnames for annGR
    seqLevels <- sub("^chr", "", seqlevels(annGR))
    seqlevels(annGR) <- seqLevels
    
    # check that all levels in GR exist in annGR, if not exclude these levels
    if (sum(!levels(seqnames(GR)) %in% levels(seqnames(annGR))) > 0) {
        TFkeepLevels <- levels(seqnames(GR)) %in% levels(seqnames(annGR))
        seqlevels(GR) <- seqlevels(GR)[TFkeepLevels]
        
    }
    
    # check that all levels in annGR exist in GR, if not exclude these levels
    if (sum(!levels(seqnames(annGR)) %in% levels(seqnames(GR))) > 0) {
        TFkeepLevels <- levels(seqnames(annGR)) %in% levels(seqnames(GR))
        seqlevels(annGR, pruning.mode="coarse") <- seqlevels(annGR)[TFkeepLevels]
    }
    
    # the seqlevels comes in different orders. This will give the correct order.
    seqlevels(GR) <- seqlevels(annGR)
    
    # add flanking regions
    lf <- flank(GR, leftFlank, start = TRUE)
    rf <- flank(GR, rightFlank, start = FALSE)
    GR <- c(lf, GR, rf)
    # find overlaps between annGR and Snps incl. flank region
    CDSInRegion <- subsetByOverlaps(annGR, GR)  #put in flankSize here when you have time ;)
    seqlengths(CDSInRegion) <- seqlengths(GR)
    CDSInRegion
    
}



#' @rdname annotation-wrappers
getCDSVector <- function(TxDb, GR, leftFlank = 0, rightFlank = 0, verbose = FALSE) {
    if (verbose) {
        cat("start CDS extraction\n")
    }
    
    CDSInRegion <- getCDSFromAnnotation(TxDb, GR, leftFlank, rightFlank)
    
    seqlevels(GR) <- seqlevels(CDSInRegion)
    
    # half-vectorized solution
    h <- findOverlaps(CDSInRegion, GR)
    CDSVec <- vector()
    for (i in 1:length(GR)) {
        CDSVec[i] <- paste(mcols(CDSInRegion[queryHits(h)[subjectHits(h) == i]])[["cds_id"]], 
            collapse = ",")
    }
    # set NAs where appropriate
    CDSVec[CDSVec == ""] <- NA
    # return list with symbols
    CDSVec
}

#' @rdname annotation-wrappers
getAnnotationDataFrame <- function(GR, strand = "+", annotationType = NULL, OrgDb = NULL, 
    TxDb = NULL, verbose = FALSE) {
    # main checks
    if (sum(!(annotationType %in% c("gene", "exon", "cds", "transcript"))) > 0) {
        stop("annotationType must be one or more of these arguments 'gene','exon','cds','transcript'")
    }
    
    if (is.null(OrgDb) & is.null(TxDb)) {
        stop("at least one of parameters OrgDb or TxDb must be used")
    }
    
    # nr of columns for return df
    ncol <- 0
    if (!is.null(OrgDb)) {
        ncol <- ncol + 1
    }
    if (!is.null(TxDb)) {
        ncol <- ncol + (length(annotationType) - 1)
    }
    
    # return dataframe
    df <- data.frame(row.names = 1:length(GR))
    
    # set strand
    strand(GR) <- strand
    
    # extract annotation
    if (!is.null(OrgDb)) {
        if ("gene" %in% annotationType) {
            gene <- getGenesVector(OrgDb = OrgDb, GR = GR, verbose = verbose)
            df[["symbol"]] <- gene
        }
        if (is.null(annotationType)) {
            gene <- getGenesVector(OrgDb = OrgDb, GR = GR, verbose = verbose)
            df[["symbol"]] <- gene
        }
    }
    if (!is.null(TxDb)) {
        if ("exon" %in% annotationType) {
            df[["exon_id"]] <- getExonsVector(TxDb = TxDb, GR = GR, verbose = verbose)
        }
        if ("transcript" %in% annotationType) {
            df[["tx_id"]] <- getTranscriptsVector(TxDb = TxDb, GR = GR, verbose = verbose)
        }
        if ("cds" %in% annotationType) {
            df[["cds_id"]] <- getCDSVector(TxDb = TxDb, GR = GR, verbose = verbose)
        }
        
        if (is.null(annotationType)) {
            df[["exon_id"]] <- getExonsVector(TxDb = TxDb, GR = GR, verbose = verbose)
            df[["tx_id"]] <- getTranscriptsVector(TxDb = TxDb, GR = GR, verbose = verbose)
            df[["cds_id"]] <- getCDSVector(TxDb = TxDb, GR = GR, verbose = verbose)
        }
    }
    df
}

#' decorateWithGenes
#' 
#' Internal function that can draw gene regions on pre-specified surfaces.
#' Necessary for the genomic-location plots.
#' 
#' The main intention of this function is to be used when plotting several bar
#' plots in the same window. This function add gene regions under the bars.
#' 
#' @param x \code{ASEset} object
#' @param genesInRegion \code{GRanges} object with gene regions. Can be
#' obtained using \code{getGenesFromAnnotation}
#' @param xlim xlim values for the pre-specified surface
#' @param ylim ylim values for the pre-specified surface
#' @param chromosome character
#' @return \code{decorateWithGenes} returns nothing, but draws genes
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{locationplot}} which is uses this
#' function internally.  \item The \code{\link{decorateWithExons}} which is
#' another similar function that \code{\link{locationplot}} uses internally.  }
#' @keywords internal
#' @examples
#' 
#'   data(ASEset)
#' 
#' @export decorateWithGenes
decorateWithGenes <- function(x, genesInRegion, xlim, ylim, chromosome) {
    
    # check the input variables
    if (class(xlim) != "integer") 
        xlim <- as.numeric(xlim)
    if (class(xlim) != "numeric") 
        stop(paste("xlim must be of class numeric, not", class(xlim)))
    if (length(xlim) != 2) 
        stop(paste("xlim must be of length 2, not", length(xlim)))
    if (class(ylim) != "integer") 
        ylim <- as.numeric(ylim)
    if (class(ylim) != "numeric") 
        stop(paste("ylim must be of class numeric, not", class(ylim)))
    if (length(ylim) != 2) 
        stop(paste("ylim must be of length 2, not", length(ylim)))
    if (class(chromosome) != "character") 
        stop(paste("chromosome must be of class character, not", class(chromosome)))
    if (length(chromosome) != 1) 
        stop(paste("chromosome must be of length 1, not", length(chromosome)))
    if (!chromosome %in% unique(seqnames(genesInRegion))) {
        if (sub("^chr", "", chromosome) %in% unique(seqnames(genesInRegion))) {
            chromosome <- sub("^chr", "", chromosome)
        } else {
            stop(paste("chromosome", chromosome, "was not found amongst the seqnames of the genesInRegion object:", 
                paste(sort(unique(seqnames(genesInRegion))), collapse = ", ")))
        }
    }
    
    
    # only work with genes on current chromosome
    genesInRegion <- genesInRegion[seqnames(genesInRegion) == chromosome]
    
    # calculate how many 'rows' have to be made available
    maxCoverage <- max(coverage(genesInRegion)[[chromosome]])
    
    # loop over all unique genes, drawing them as specified
    uniqueGenes <- unique(mcols(genesInRegion)[["Symbol"]])
    
    for (i in 1:length(uniqueGenes)) {
        # getting the name of the gene and the height on the Y-axis
        genesymbol <- uniqueGenes[i]
        yPos <- ylim[1] + (i - 1)%%maxCoverage * ((ylim[2] - ylim[1])/maxCoverage)
        
        # this block checks for double instances and just arbitrarily take the first
        # (typically miRNAs with two locations)
        if (sum(mcols(genesInRegion)[["Symbol"]] %in% genesymbol) > 1) {
            geneData <- genesInRegion[which(mcols(genesInRegion)[["Symbol"]] %in% 
                genesymbol)[1], ]
        } else {
            geneData <- genesInRegion[which(mcols(genesInRegion)[["Symbol"]] %in% 
                genesymbol), ]
        }
        
        # draw and label
        start <- max(c(xlim[1] - (xlim[2] - xlim[1])/10, start(geneData)))
        end <- min(c(xlim[2] + (xlim[2] - xlim[1])/10, end(geneData)))
        lines(x = c(start, end), y = c(yPos, yPos), lwd = 2)
        text(x = start + (end - start)/2, y = yPos + (ylim[2] - ylim[1])/6, label = genesymbol, 
            cex = 0.8)
    }
}



#' decorateWithExons
#' 
#' Internal function that can draw gene regions on pre-specified surfaces.
#' Necessary for the genomic-location plots.
#' 
#' The main intention of this function is to be used when plotting several bar
#' plots in the same window. This function add gene regions under the bars.
#' 
#' @param x \code{ASEset} object
#' @param exonsInRegion \code{GRanges} object with generegions. Can be obtained
#' using \code{getExonsFromAnnotation}. Must contain a column 'tx_name'
#' @param xlim xlim values for the pre-specified surface
#' @param ylim ylim values for the pre-specified surface
#' @param chromosome character
#' @return \code{decorateWithExons} returns nothing, but draws genes
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{locationplot}} which is uses this
#' function internally.  \item The \code{\link{decorateWithGenes}} which is
#' another similar function that \code{\link{locationplot}} uses internally.  }
#' @keywords internal
#' @examples
#' 
#'   data(ASEset)
#' 
#' @export decorateWithExons
decorateWithExons <- function(x, exonsInRegion, xlim, ylim, chromosome) {
    
    # check the input variables
    if (class(exonsInRegion) != "GRanges") 
        stop(paste("exonsInRegion must be of class GRanges, not", class(exonsInRegion)))
    if (!"tx_name" %in% colnames(mcols(exonsInRegion))) 
        stop("exonsInRegion must contain an mcol variable named 'tx_name'")
    if (class(xlim) != "integer") 
        xlim <- as.numeric(xlim)
    if (class(xlim) != "numeric") 
        stop(paste("xlim must be of class numeric, not", class(xlim)))
    if (length(xlim) != 2) 
        stop(paste("xlim must be of length 2, not", length(xlim)))
    if (class(ylim) != "integer") 
        ylim <- as.numeric(ylim)
    if (class(ylim) != "numeric") 
        stop(paste("ylim must be of class numeric, not", class(ylim)))
    if (length(ylim) != 2) 
        stop(paste("ylim must be of length 2, not", length(ylim)))
    if (class(chromosome) != "character") 
        stop(paste("chromosome must be of class character, not", class(chromosome)))
    if (length(chromosome) != 1) 
        stop(paste("chromosome must be of length 1, not", length(chromosome)))
    if (!chromosome %in% unique(seqnames(exonsInRegion))) {
        if (sub("^chr", "", chromosome) %in% unique(seqnames(exonsInRegion))) {
            chromosome <- sub("^chr", "", chromosome)
        } else {
            stop(paste("chromosome", chromosome, "was not found amongst the seqnames of the exonInRegion object:", 
                paste(sort(unique(seqnames(exonsInRegion))), collapse = ", ")))
        }
    }
    # only work with exons on current chromosome
    exonsInRegion <- exonsInRegion[seqnames(exonsInRegion) == chromosome]
    
    
    # calculate how many 'rows' have to be made available (corresponding to the
    # number of unique transcripts in exonsInRegion
    uniqueGenes <- unique(unlist(as.list(mcols(exonsInRegion)[["tx_name"]])))
    maxCoverage <- length(uniqueGenes)
    
    
    for (i in 1:length(uniqueGenes)) {
        # getting the name of the gene and the height on the Y-axis
        tx_name <- uniqueGenes[i]
        yPos <- ylim[1] + (i - 1)%%maxCoverage * ((ylim[2] - ylim[1])/maxCoverage)
        
        # extracting and drawing all exons for each transcript
        exonsInTranscript <- which(sapply(as.list(mcols(exonsInRegion)[["tx_name"]]), 
            function(x) {
                tx_name %in% x
            }))
        for (exonInTranscript in exonsInTranscript) {
            x1 <- start(exonsInRegion[exonInTranscript])
            x2 <- end(exonsInRegion[exonInTranscript])
            lines(x = c(x1, x2), y = c(yPos, yPos), lwd = 2)
        }
        
        
        # draw a thin connecting line for each transcript
        end <- max(end(exonsInRegion[exonsInTranscript]))
        start <- min(start(exonsInRegion[exonsInTranscript]))
        lines(x = c(start, end), y = c(yPos, yPos), lwd = 0.5)
        
        # label with the tx_name
        start <- max(c(start, xlim[1]))
        end <- min(c(end, xlim[2]))
        text(x = start + (end - start)/2, y = yPos + (ylim[2] - ylim[1])/6, label = tx_name, 
            cex = 0.8)
        
    }
} 
