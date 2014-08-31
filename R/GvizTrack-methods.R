#'@include locationplot-methods.R
NULL

#' ASEset-gviztrack ASEset objects
#' 
#' plotting ASE effects over a specific genomic region
#' 
#' For information of how to use these tracks in more ways, visit the Gviz
#' package manual.
#' 
#' @name ASEset-gviztrack
#' @rdname ASEset-gviztrack
#' @aliases ASEset-gviztrack CoverageDataTrack ASEDAnnotationTrack
#' CoverageDataTrack,ASEset-method ASEDAnnotationTrack,ASEset-method
#' @docType methods
#' @param x an ASEset object.
#' @param GR genomic range of plotting
#' @param type 'fraction' or 'count'
#' @param strand '+','-'. This argument determines which strand is plotted.
#' @param mainVec vector of text for the main of each plot
#' @param BamList GAlignmnentsList object of reads from the same genomic region
#' as the ASEset
#' @param start start position of reads to be plotted
#' @param trackName name of track (ASEDAnnotationTrack)
#' @param trackNameVec names of tracks (CoverageDataTrack)
#' @param end end position of reads to be plotted
#' @param gpar graphical parameters for each small plot
#' @param verbose Setting \code{verbose=TRUE} gives details of procedure during
#' function run
#' @param ... arguments passed on to barplot function
#' @author Jesper R. Gadin
#' @seealso \itemize{ \item The \code{\link{ASEset}} class which the functions
#' can be called up on.}
#' @keywords ASEDAnnotationTrack CoverageDataTrack
#' @examples
#' 
#' data(ASEset)
#' x <- ASEset[,1:2]
#' r <- reads[1:2]
#' genome(x) <- 'hg19'
#' seqlevels(r) <- seqlevels(x)
#' 
#' GR <- GRanges(seqnames=seqlevels(x),ranges=IRanges(start=min(start(x)),end=max(end(x))),strand='+', genome=genome(x))
#' 
#' deTrack <- ASEDAnnotationTrack(x, GR=GR, type='fraction',strand='+')
#' covTracks <- CoverageDataTrack(x,BamList=r,strand='+') 
#' 
#' lst <- c(deTrack,covTracks)
#' 
#' sizes <- c(0.5,rep(0.5/length(covTracks),length(covTracks)))
#' #temporarily do not run this function 
#' #plotTracks(lst, from=min(start(x)), to=max(end(x)), 
#' #sizes=sizes, col.line = NULL, showId = FALSE, main='mainText', 
#' #cex.main=1, title.width=1, type='histogram')
#' 
#'
#' @importFrom Gviz DetailsAnnotationTrack
#' @importFrom Gviz AnnotationTrack
#' @importFrom Gviz DataTrack
#' @importFrom Gviz plotTracks
#' 
#' @importClassesFrom Gviz DataTrack
#' @importClassesFrom Gviz AnnotationTrack
#' @importClassesFrom Gviz DetailsAnnotationTrack
#'
#' @exportMethod ASEDAnnotationTrack
#' @exportMethod CoverageDataTrack
NULL

# @export plotTracks

##' @rdname ASEset-gviztrack
# setGeneric('plotTracks')

#' @rdname ASEset-gviztrack
setGeneric("ASEDAnnotationTrack", function(x, GR = rowData(x), type = "fraction", 
    strand = "+", mainVec = vector(), trackName = paste("deTrack", type), verbose = TRUE, 
    gpar = list(), ...) {
    standardGeneric("ASEDAnnotationTrack")
})

setMethod("ASEDAnnotationTrack", signature(x = "ASEset"), function(x, GR = rowData(x), 
    type = "fraction", strand = "+", mainVec = rep("", nrow(x)), trackName = paste("deTrack", 
        type), verbose = TRUE, gpar = list(), ...) {
    
    # change to '*' if(strand=='nonStranded'){strand <- '*'} not possile as long as
    # nonStranded is an option
    
    # check genome
    if (is.null(genome(x)) | is.na(genome(x))) {
        stop(paste("genome have be set for object x", "e.g. genome(x) <- \"hg19\" "))
    }
    
    # check seqnames has length=1
    if (!(length(seqlevels(x)) == 1)) {
        stop("This function can only use objects with one seqlevel")
    }
    
    if (sum(strand == "+" | strand == "-") == 0) {
        stop("strand must be plus or minus at the moment")
    }
    if (!nrow(x) == 1) {
        GR <- GRanges(seqnames = seqlevels(x), ranges = IRanges(start = min(start(x)), 
            end = max(end(x))), strand = strand, genome = genome(x))
        # if(sum(width(reduce(GR)))==1 ){ GR <- flank(GR,2,both=TRUE) }
    }
    
    # check gpar, set default if not set by user
    gparDefault <- list(ylab = "", xlab = "", deAnnoPlot = TRUE)
    gparSet <- names(gparDefault) %in% names(gpar)
    gpar <- mapply(gparDefault, gparSet, names(gparDefault), FUN = function(x, y, 
        z, gpar) {
        if (y) {
            gpar[[z]]
        } else {
            x
        }
        
    }, MoreArgs = list(gpar = gpar))
    
    ranges <- rowData(x)
    
    colnames(x) <- 1:ncol(x)
    
    details <- function(identifier, ...) {
        
        # an option to use more print(list(...))
        
        type <- get("type", envir = AllelicImbalance.extra)
        arank <- get("arank", envir = AllelicImbalance.extra)
        afraction <- get("afraction", envir = AllelicImbalance.extra)
        acounts <- get("acounts", envir = AllelicImbalance.extra)
        amainVec <- get("amainVec", envir = AllelicImbalance.extra)
        # gparam <- get('gpar',envir=AllelicImbalance.extra)
        
        
        if (type == "fraction") {
            print(barplotLatticeFraction(identifier, afraction, arank, amainVec, 
                ...), newpage = FALSE, prefix = "plot")
            
        } else if (type == "count") {
            print(barplotLatticeCounts(identifier, acounts, arank, amainVec, ...), 
                newpage = FALSE, prefix = "plot")
        }
        
    }
    
    # pick out plot data from ASEset strand='+'
    AllelicImbalance.extra <- new.env(parent = emptyenv())
    amainVec <- mainVec
    
    assign("acounts", alleleCounts(x, strand = strand), envir = AllelicImbalance.extra)
    assign("arank", arank(x, strand = strand), envir = AllelicImbalance.extra)
    assign("afraction", fraction(x, strand = strand), envir = AllelicImbalance.extra)
    assign("type", type, envir = AllelicImbalance.extra)
    assign("amainVec", amainVec, envir = AllelicImbalance.extra)
    
    
    # plot the fraction
    deTrack <- AnnotationTrack(range = ranges, genome = genome(x), id = rownames(x), 
        name = trackName, stacking = "squish", fun = details, detailsFunArgs = gpar)
    deTrack
})

#' @rdname ASEset-gviztrack
setGeneric("CoverageDataTrack", function(x, GR = rowData(x), BamList = NULL, strand = NULL, 
    start = NULL, end = NULL, trackNameVec = NULL, verbose = TRUE, ...) {
    standardGeneric("CoverageDataTrack")
})



setMethod("CoverageDataTrack", signature(x = "ASEset"), function(x, GR = NULL, BamList = NULL, 
    strand = NULL, start = NULL, end = NULL, trackNameVec = NULL, verbose = TRUE, 
    ...) {
    # GR is not in use atm. Missing is a subset of the return matrix based on the GR
    # values.
    
    if (!is.null(strand)) {
        if (strand == "+") {
            pstrand = TRUE
        } else if (strand == "-") {
            mstrand = TRUE
        } else {
            stop("strand has to be '+' or '-' if not NULL\n")
        }
    } else {
        stop("strand has to be '+' or '-' if not NULL\n")
    }
    
    # check genome
    if (is.null(genome(x)) | is.na(genome(x))) {
        stop(paste("genome have be set for object x", "e.g. genome(x) <- \"hg19\" "))
    }
    
    # check seqnames has length=0
    if (!(length(seqlevels(x)) == 1)) {
        stop("This function can only use objects with one seqlevel")
    }
    
    if (!is.null(strand)) {
        if (strand == "+") {
            pstrand = TRUE
        } else if (strand == "-") {
            mstrand = TRUE
        } else {
            stop("strand has to be '+' or '-' if not NULL\n")
        }
    }
    
    GR <- GRanges(seqnames = seqlevels(x), ranges = IRanges(start = min(start(x)), 
        end = max(end(x))), strand = strand)
    
    # start <- start(GR) end <- end(GR) chr <- seqnames(GR)
    
    if (is.null(BamList)) {
        stop("must include GappedAlignmentsList as BamList ")
    } else {
        # check that only one chromosome is present
        if (!length(seqlevels(BamList)) == 1) {
            stop("can only be one seq level\n")
        }
        covMatList <- coverageMatrixListFromGAL(BamList, strand)
    }
    
    trackList <- list()
    
    mat <- covMatList[["mat"]]
    start <- covMatList[["start"]]
    end <- covMatList[["end"]]
    
    if (is.null(trackNameVec)) {
        trackNameVec[1:ncol(x)] <- colnames(x)
    } else {
        if (!length(trackNameVec) == nrow(x)) {
            stop("length of trackNameVec must be equal to cols in (ASEset)")
        }
    }
    
    # prepare Gviz dtracks
    for (j in 1:nrow(mat)) {
        trackList[[length(trackList) + 1]] <- DataTrack(data = mat[j, ], start = start:end, 
            width = 1, chromosome = seqlevels(x), genome = genome(x), name = trackNameVec[j], 
            type = "s")
    }
    
    trackList
}) 
