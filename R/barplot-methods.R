
#' barplot ASEset objects
#' 
#' Generates barplots for ASEset objects. Two levels of plotting detail are
#' provided: a detailed barplot of read counts by allele useful for fewer
#' samples and SNPs, and a less detailed barplot of the fraction of imbalance,
#' useful for more samples and SNPs.
#' 
#' \code{filter.pValue.fraction} is intended to remove p-value annotation with
#' very large difference in frequency, which could just be a sequencing
#' mistake. This is to avoid p-values like 1e-235 or similar.
#' 
#' \code{sampleColour}User specified colours, either given as named colours
#' ('red', 'blue', etc) or as hexadecimal code. Can be either length 1 for all
#' samples, or else of a length corresponding to the number of samples for
#' individual colouring.
#' 
#' @name ASEset-barplot
#' @rdname ASEset-barplot
#' @aliases ASEset-barplot barplot barplot,ASEset-method
#' @docType methods
#' @param height An \code{ASEset} object
#' @param type 'count' or 'fraction'
#' @param sampleColour User specified colours
#' @param legend Display legend
#' @param pValue Display p-value
#' @param strand four options, '+', '-', 'both' or '*'
#' @param testValue if set, a matrix or vector with user p-values
#' @param testValue2 if set, a matrix or vector with user p-values
#' @param OrgDb an OrgDb object which provides annotation
#' @param TxDb a TxDb object which provides annotation
#' @param annotationType select one or more from
#' 'gene','exon','transcript','cds'.
#' @param main text to use as main label
#' @param ylim set plot y-axis limit
#' @param yaxis wheter the y-axis is to be displayed or not
#' @param xaxis wheter the x-axis is to be displayed or not
#' @param ylab showing labels for the tic marks
#' @param xlab showing labels for the tic marks
#' @param legend.colnames gives colnames to the legend matrix
#' @param las.ylab orientation of ylab text
#' @param las.xlab orientation of xlab text
#' @param cex.main set main label size (max 2)
#' @param cex.pValue set pValue label size
#' @param cex.ylab set ylab label size
#' @param cex.xlab set xlab label size
#' @param cex.legend set legend label size
#' @param cex.annotation size of annotation text
#' @param ypos.annotation relative ypos for annotation text
#' @param annotation.interspace space between annotation text
#' @param legend.interspace set legend space between fills and text
#' @param add \code{boolean} indicates if a new device should be started
#' @param lowerLeftCorner integer that is only useful when \code{add}=TRUE
#' @param size Used internally by locationplot. Rescales each small barplot
#' window
#' @param addHorizontalLine adds a horizontal line that marks the default
#' fraction of 0.5 - 0.5
#' @param add.frame \code{boolean} to give the new plot a frame or not
#' @param filter.pValue.fraction \code{numeric} between 0 and 1 that filter
#' away pValues where the main allele has this frequency.
#' @param legend.fill.size size of the fill/boxes in the legend (default:NULL)
#' @param verbose Makes function more talkative
#' @param top.allele.criteria 'maxcount', 'ref' or 'phase'
#' @param ... for simpler generics when extending function
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{ASEset}} class which the barplot
#' function can be called up on.  }
#' @keywords barplot
#' @examples
#' 
#' data(ASEset)
#' barplot(ASEset[1])
#'
#' @importFrom graphics plot
#' @importFrom graphics barplot
#' @exportMethod barplot
NULL

#' @rdname ASEset-barplot
setGeneric("barplot")

#' @rdname ASEset-barplot
setMethod("barplot", signature(height = "ASEset"), function(height, type = "count", 
    sampleColour = NULL, legend = TRUE, pValue = TRUE, strand = "*", testValue = NULL, 
    testValue2 = NULL, OrgDb = NULL, TxDb = NULL, annotationType = c("gene", "exon", 
        "transcript"), main = NULL, ylim = NULL, yaxis = TRUE, xaxis = FALSE, ylab = TRUE, 
    xlab = TRUE, legend.colnames = "", las.ylab = 1, las.xlab = 2, cex.main = 1, 
    cex.pValue = 0.7, cex.ylab = 0.7, cex.xlab = 0.7, cex.legend = 0.6, add = FALSE, 
    lowerLeftCorner = c(0, 0), size = c(1, 1), addHorizontalLine = 0.5, add.frame = TRUE, 
    filter.pValue.fraction = 0.99,  legend.fill.size=1, legend.interspace=1, verbose = FALSE, 
	top.allele.criteria="maxcount", cex.annotation=0.7, ypos.annotation=0, annotation.interspace=1,
	...) {
    
    # catch useful graphical parameters that can be used to later add onto plot. This
    # list will be retireved by using 'glst <- barplot(x)'
    graphParamList <- list()
    
    x <- height
    
    # check type
    okPlotTypes <- c("fraction", "count")
    if (!type %in% okPlotTypes) 
        stop(paste("type can't be '", type, "' - it should be one of these: ", paste(okPlotTypes, 
            collapse = ", "), sep = ""))
    
    # check verbose
    if (class(verbose) != "logical") 
        stop(paste("verbose should be of class logical, not", class(verbose)))
    if (length(verbose) != 1) 
        stop(paste("verbose should be of length 1, not", length(verbose)))
    
    
    # check strand
    okStrandTypes <- c("both", "+", "-", "*")
    if (!strand %in% okStrandTypes) 
        stop(paste("strand can't be '", strand, "' - it should be one of these: ", 
            paste(okStrandTypes, collapse = ", "), sep = ""))
    
    # check if strand type is present in object
    if (strand == "+") {
        el <- "countsPlus"
        if (!(el %in% names(assays(x)))) {
            stop("strand is not present as assay in ASEset object")
        }
    } else if (strand == "-") {
        el <- "countsMinus"
        if (!(el %in% names(assays(x)))) {
            stop("strand is not present as assay in ASEset object")
        }
    } else if (strand == "*") {
        el <- "countsUnknown"
        if (!(el %in% names(assays(x)))) {
            stop("strand is not present as assay in ASEset object")
        }
    } else if (strand == "both") {
        el <- "countsPlus"
        if (!(el %in% names(assays(x)))) {
            stop("strand is not present as assay in ASEset object")
        }
        el <- "countsMinus"
        if (!(el %in% names(assays(x)))) {
            stop("strand is not present as assay in ASEset object")
        }
    } else {
        stop("unknown strand option")
    }
    
    # check legend
    if (!is.logical(legend)) {
        stop("legend has to be logical, TRUE or FALSE")
    }
    if (!length(legend) == 1) {
        stop("legend has to be of length 1")
    }
    
    # check pValue
    if (!is.logical(pValue)) {
        stop("pValue has to be logical, TRUE or FALSE")
    }
    if (!length(pValue) == 1) {
        stop("pValue has to be of length 1")
    }
    
    
    # check ylim
    if (!is.null(ylim)) {
    }
    
    # check testValue
    if (!is.null(testValue)) {
        if (!class(testValue) == "numeric" & !class(testValue) == "matrix") {
            stop("wrong class for testValue")
        }
        if (!class(testValue) == "matrix") {
            if (!sum(dim(testValue) == rev(dim(x))) == 2) {
                stop("wrong dimensions of testValue")
            }
        }
        if (!class(testValue) == "numeric") {
            if (!length(testValue) == ncol(x)) {
                stop("wrong dimensions of testValue")
            }
        }
    }
    # check testValue2
    if (!is.null(testValue2)) {
        if (!class(testValue2) == "numeric" & !class(testValue2) == "matrix") {
            stop("wrong class for testValue")
        }
        if (!class(testValue2) == "matrix") {
            if (!sum(dim(testValue2) == rev(dim(x))) == 2) {
                stop("wrong dimensions of testValue")
            }
        }
        if (!class(testValue2) == "numeric") {
            if (!length(testValue2) == ncol(x)) {
                stop("wrong dimensions of testValue")
            }
        }
    }
    
    # check testValue combinations
    if (!strand == "both" & !is.null(testValue2)) {
        stop("testValue2 is supposed to only be \n\t\t\tspecified when using strand='both', and then be the values for the minus strand.")
    }
    
    # check ylim
    if (!is.null(ylim)) {
        if (type == "fraction") {
            warning("an ylim argument given to barplot with type fraction is ignored")
            ylim <- NULL
        } else {
            if (class(ylim) != "numeric") 
                stop(paste("if given, ylim must be of class numeric, not", class(ylim)))
            if (length(ylim) != 2) 
                stop(paste("if given, ylim must be of length 2, not", length(ylim)))
            if (ylim[2] < ylim[1]) 
                stop(paste("if ylim is given, ylim[2] must be larger than ylim[1]"))
        }
    }
    
    # set easy accessor to samples
    samples <- colnames(x)
    # set easy accessor to snps
    snps <- rownames(x)
    
    
    # check sampleColour
    if (is.null(sampleColour)) {
        sampleColour <- rep("dodgerblue", length(samples))
    } else {
        if (class(sampleColour) != "character") 
            stop(paste("if given, sampleColour should be of class character, not", 
                class(sampleColour)))
        if (length(sampleColour) == 1) {
            sampleColour <- rep(sampleColour, length(samples))
        } else {
            if (length(sampleColour) != length(samples)) 
                stop(paste("if given, sampleColour should be of length 1 or length", 
                  length(samples), "(corresponding to the number of samples)"))
        }
        hexadecimals <- sampleColour[grep("^#[0-9A-F]{6}", sampleColour)]
        missing <- sampleColour[!(sampleColour %in% colors() | sampleColour %in% 
            hexadecimals)]
        if (length(missing) > 0) {
            stop(paste("The following", length(unique(missing)), "colour(s) from sampleColour were not recognised:", 
                paste(sort(unique(missing)), collapse = ", ")))
        }
    }
    
    # check OrgDb
    if (!is.null(OrgDb)) {
        if (!class(OrgDb) == "OrgDb") {
            stop("class of OrgDb has to be an OrgDb class")
        }
    }
    if (!is.null(TxDb)) {
        if (!class(TxDb) == "TxDb") {
            stop("class of TxDb has to be a TxDb class")
        }
    }
    
    # check annotationType
    okAnnotationTypes <- c("gene", "exon", "cds", "transcript")
    if (!sum(annotationType %in% okAnnotationTypes) == length(annotationType)) {
        stop(paste("type can't be '", annotationType, "' - it should be one of these: ", 
            paste(okAnnotationTypes, collapse = ", "), sep = ""))
    }
    
    # check main
    if (!is.null(main)) {
        if (!class(main) == "character") {
            stop("main has to be of class character")
        }
    }
    
    # check yaxis
    if (!is.logical(yaxis)) {
        stop("yaxis has to be logical, TRUE or FALSE")
    }
    if (!length(yaxis) == 1) {
        stop("yaxis has to be of length 1")
    }
    
    # check xaxis
    if (!is.logical(xaxis)) {
        stop("xaxis has to be logical, TRUE or FALSE")
    }
    if (!length(xaxis) == 1) {
        stop("xaxis has to be of length 1")
    }
    
    # check ylab
    if (!is.logical(ylab)) {
        stop("ylab has to be logical, TRUE or FALSE")
    }
    if (!length(ylab) == 1) {
        stop("ylab has to be of length 1")
    }
    
    # check xlab
    if (!is.logical(xlab)) {
        stop("xlab has to be logical, TRUE or FALSE")
    }
    if (!length(xlab) == 1) {
        stop("xlab has to be of length 1")
    }
    
    # check legend.colnames
    if (!class(legend.colnames) == "character") {
        stop("legend.colnames has to be of class character")
    }
    
    # check las.ylab
    if (!class(las.ylab) == "numeric" & !class(las.ylab) == "integer") {
        stop("las.ylab has to be of class numeric or integer")
    }
    if (!las.ylab == 1 & !las.ylab == 2) {
        stop("las.ylab has to be a number, either 1 or 2")
    }
    
    # check las.xlab
    if (!class(las.xlab) == "numeric" & !class(las.xlab) == "integer") {
        stop("las.xlab has to be of class numeric or integer")
    }
    if (!las.xlab == 1 & !las.xlab == 2) {
        stop("las.xlab has to be a number, either 1 or 2")
    }
    
    # check cex all of them
    if (!class(cex.main) == "numeric") {
        if (cex.main <= 0 | !length(cex.main) == 1) 
            stop("cex.main has to be of class numeric, and have length 1, and a value above 0 ")
    }
    if (!class(cex.pValue) == "numeric") {
        if (cex.pValue > 1 | cex.pValue <= 0 | !length(cex.pValue) == 1) 
            stop("cex.pValue has to be of class numeric, and have length 1, and a value above 0")
    }
    if (!class(cex.ylab) == "numeric") {
        if (cex.ylab <= 0 | !length(cex.ylab) == 1) 
            stop("cex.ylab has to be of class numeric, and have length 1, and a value above 0 ")
    }
    if (!class(cex.xlab) == "numeric") {
        if (cex.xlab <= 0 | !length(cex.xlab) == 1) 
            stop("cex.xlab has to be of class numeric, and have length 1, and a value above 0 ")
    }
    if (!class(cex.legend) == "numeric") {
        if (cex.legend <= 0 | !length(cex.legend) == 1) 
            stop("cex.legend has to be of class numeric, and have length 1, and a value above 0 ")
    }
    
    # check add
    if (!is.logical(add)) {
        stop("add has to be logical, TRUE or FALSE")
    }
    if (!length(add) == 1) {
        stop("add has to be of length 1")
    }
    
    # check lowerLeftCorner
    if (!class(lowerLeftCorner) == "numeric" & !class(lowerLeftCorner) == "integer") {
        stop("lowerLeftCorner has to be of class numeric or integer")
    }
    if (!length(lowerLeftCorner) == 2) {
        stop("lowerLeftCorner has be of length 2")
    }
    
    # check size
    if (!class(size) == "numeric" & !class(size) == "integer") {
        stop("size has to be of class numeric or integer")
    }
    if (!length(size) == 2) {
        stop("size has be of length 2")
    }
    
    # check addHorizontalLine
    if (!is.null(addHorizontalLine)) {
        if (!class(addHorizontalLine) == "numeric") {
            stop("addHorizontalLine has to be either NULL or of class numeric")
        }
        if (!length(addHorizontalLine) == 1) {
            stop("addHorizontalLine has be of length 1")
        }
        if (addHorizontalLine > 1 | addHorizontalLine <= 0) {
            stop("addHorizontalLine has be numeric between 0 and 1")
        }
    }
    
    # check add.frame
    if (!is.logical(add.frame)) {
        stop("add.frame has to be logical, TRUE or FALSE")
    }
    if (!length(add.frame) == 1) {
        stop("add.frame has to be of length 1")
    }
    
    # check filter.pValue.fraction=0.99
    if (!class(filter.pValue.fraction) == "numeric") {
        stop("filter.pValue.fraction has to be of class numeric")
    }
    if (!length(filter.pValue.fraction) == 1) {
        stop("filter.pValue.fraction has to be of length 1")
    }
    if (filter.pValue.fraction < 0 | filter.pValue.fraction > 1) {
        stop("filter.pValue.fraction has to be a value between 0 and 1")
    }
    
    # make annotation dataframes
    if ((!is.null(OrgDb) | !is.null(TxDb)) & (strand == "+" |  
        strand == "*" | strand == "both")) {
        annDfPlus <- getAnnotationDataFrame(rowRanges(x), strand = "+", annotationType = annotationType, 
            OrgDb = OrgDb, TxDb = TxDb, verbose = verbose)
    } else {
        annDfPlus <- NULL
    }
    if ((!is.null(OrgDb) | !is.null(TxDb)) & (strand == "-" |  
        strand == "*" | strand == "both")) {
        annDfMinus <- getAnnotationDataFrame(rowRanges(x), strand = "-", annotationType = annotationType, 
            OrgDb = OrgDb, TxDb = TxDb, verbose = verbose)
    } else {
        annDfMinus <- NULL
    }
    
    # creating a background from foreground colour - basically just turning down the
    # saturation a little compared to the foreground color except for very
    # low-saturated fg's which are made more saturated, and very dark fg's which are
    # made lighter
    hsvColoursBg <- hsvColoursFg <- rgb2hsv(col2rgb(sampleColour))
    veryLowSaturation <- hsvColoursFg["s", ] < 0.3
    veryDark <- hsvColoursFg["v", ] < 0.2
    remaining <- !veryLowSaturation & !veryDark
    hsvColoursBg["s", veryLowSaturation] <- hsvColoursFg["s", veryLowSaturation] + 
        0.3
    hsvColoursBg["v", veryDark] <- hsvColoursFg["v", veryDark] + 0.2
    hsvColoursBg["s", remaining] <- hsvColoursFg["s", remaining] - 0.3
    fgCol <- hsv(hsvColoursFg["h", ], hsvColoursFg["s", ], hsvColoursFg["v", ])
    bgCol <- hsv(hsvColoursBg["h", ], hsvColoursBg["s", ], hsvColoursBg["v", ])
    
    
    if (type == "count") {
        # the intention of this is to be able to create barplots which show both the
        # amount of counts and the distribution between alleles this is the least
        # digested format and consequently produces the most cluttered plots. It is
        # therefore intended for subsets, few snps, few samples etc in these smaller
        # subsets it is useful to get information on higher expression values at certain
        # points of genes etc.
        
        if (strand == "+" | strand == "-" | strand == "*" ) {
            for (i in 1:length(snps)) {
                # getting revelant data
                snp <- snps[i]
                tmp <- t(alleleCounts(x, strand)[[snp]][samples, , drop = FALSE])
                
                if ((ncol(tmp) * nrow(tmp)) == sum(tmp == 0)) {
                  zeroRows <- TRUE
                } else {
                  zeroRows <- FALSE
                }
                
                # alleles are
                n <- names(sort(apply(tmp, 1, sum), decreasing = TRUE)[1:2])
                # theor total count
                tc <- sort(apply(tmp, 1, sum), decreasing = TRUE)[1:2]
                
                # extracting only major and minor allele in correctly transposed format and
                # ordering colours accordingly
                height <- tmp[match(n, rownames(tmp)), , drop = FALSE]
                # change all zero rows to X names
                rownames(height)[apply(height, 1, sum) == 0] <- "X"
                
                # warn if the two remaining alleles have too many counts
                majorAndMinorFraction <- sum(height)/sum(tmp)
                if (verbose & majorAndMinorFraction < 0.9) {
                  message(paste("Snp", snp, "was possible tri-allelic, but only two most frequent alleles were plotted. Counts:\n", 
                    paste(paste(names(apply(tmp, 1, sum)), apply(tmp, 1, sum), sep = "="), 
                      collapse = ", ")))
                }
                
                
                ## make bars section
                extraspaceOnTop <- 1.15
                if (is.null(ylim)) {
                  if (!zeroRows) {
                    ylimHere <- c(0, max(abs(range(height, na.rm = TRUE))) * extraspaceOnTop)
                  } else {
                    ylimHere <- c(0, 10)
                  }
                } else {
                  ylimHere <- ylim
                }
                
                # init frame
                if (!add) {
                  size <- c(1, 1)
                  lowerLeftCorner <- c(0, 0)
                  plot.default(NULL, xlim = c(0, 1), ylim = c(0, 1), ylab = "", xlab = "", 
                    xaxt = "n", yaxt = "n", frame.plot = FALSE)
                }
                
                
                # set so the ylimHere is at least 10.
                if (ylimHere[2] < 10) {
                  ylimHere[2] <- 10
                }
                
                # make tics and labels prettier
                ylabels <- pretty(ylimHere[1]:ylimHere[2])
                if (strand == "both") {
                  ylimHere[1] <- -max(abs(ylabels))
                  ylimHere[2] <- max(abs(ylabels))
                } else {
                  ylimHere[1] <- min(ylabels)
                  ylimHere[2] <- max(ylabels)
                }
                
                # all x positions
                xPointsAll <- seq(lowerLeftCorner[1], lowerLeftCorner[1] + size[1], 
                  length.out = length(height) * 1.5 + 2)[2:(length(height) * 1.5 + 
                  1)]
                widthsOfRectangles <- rep(size[1]/(length(height) * 1.5 + 1), length(height)/2)
                # normalize toPlot-values to a range of c(0,size[2])
                heightsOfRectangles <- ((height - ylimHere[1])/(ylimHere[2] - ylimHere[1]) * 
                  size[2])
                
                
                # take every third element (start at 2)
                xPoints1 <- xPointsAll[seq(1, length(height) * 1.5, by = 3)]
                yPoints <- lowerLeftCorner[2] + heightsOfRectangles[1, ]/2
                
                # alleleColors1 <- rep(fgCol[1],9)
                alleleColors1 <- fgCol
                symbols(x = xPoints1, y = yPoints, bg = alleleColors1, rectangles = matrix(ncol = 2, 
                  c(widthsOfRectangles, heightsOfRectangles[1, ])), add = TRUE, inches = FALSE)
                
                # take every third element (start at 3)
                xPoints2 <- xPointsAll[seq(2, length(height) * 1.5, by = 3)]
                yPoints <- lowerLeftCorner[2] + heightsOfRectangles[2, ]/2
                
                alleleColors2 <- rep(bgCol[1], 9)
                alleleColors2 <- bgCol
                widthsOfRectangles <- rep(size[1]/(length(height) * 1.5 + 1), length(height)/2)
                symbols(x = xPoints2, y = yPoints, bg = alleleColors2, rectangles = matrix(ncol = 2, 
                  c(widthsOfRectangles, heightsOfRectangles[2, ])), add = TRUE, inches = FALSE)
                
                ## make axis and labels section
                
                # y-axis
                yat <- (ylabels + abs(min(ylabels)))/(sum(abs(min(ylabels)) + max(ylabels)))
                graphParamList[["yat"]] <- yat
                if (yaxis) {
                  # vertical line
                  lines(c(lowerLeftCorner[1], lowerLeftCorner[1]), c(lowerLeftCorner[2], 
                    lowerLeftCorner[2] + 1))
                  # horizontal lines
                  invisible(mapply(FUN = function(x, y) {
                    lines(c(x, x - size[1] * 0.01), c(y, y))
                  }, rep(lowerLeftCorner[1], length(yat)), (lowerLeftCorner[2] + 
                    yat)))
                }
                if (ylab) {
                  if (las.ylab == 2) {
                    text(lowerLeftCorner[1] - size[1] * 0.02, yat, labels = abs(ylabels), 
                      xpd = TRUE, srt = 90, adj = c(1, 0.5), cex = cex.ylab)
                  }
                  if (las.ylab == 1) {
                    text(lowerLeftCorner[1] - size[1] * 0.02, yat, labels = abs(ylabels), 
                      xpd = TRUE, srt = 0, adj = c(1, 0.5), cex = cex.ylab)
                  }
                }
                # x-axis labels
                xat <- xPoints1 + (xPoints2 - xPoints1)/2
                graphParamList[["xat"]] <- xat
                
                if (xaxis) {
                  # vertical line
                  lines(c(lowerLeftCorner[1], lowerLeftCorner[1] + 1), c(lowerLeftCorner[2], 
                    lowerLeftCorner[2]))
                  # horizontal lines
                  invisible(mapply(FUN = function(x, y) {
                    lines(c(x, x), c(y, y - size[1] * 0.01))
                  }, (lowerLeftCorner[2] + xat), rep(lowerLeftCorner[1], length(xat))))
                }
                
                if (xlab) {
                  if (las.xlab == 2) {
                    text(xat, lowerLeftCorner[2] - size[2] * 0.02, labels = samples, 
                      xpd = TRUE, srt = 90, adj = c(1, 0.5), cex = cex.xlab)
                  }
                  if (las.xlab == 1) {
                    text(xat, lowerLeftCorner[2] - size[2] * 0.02, labels = samples, 
                      xpd = TRUE, srt = 0, adj = c(1, 0.5), cex = cex.xlab)
                  }
                }
                
                if (is.null(main)) {
                  text(x = lowerLeftCorner[1] + size[1] * 0.5, y = lowerLeftCorner[2] + 
                    size[2] * 1.05, labels = paste(snp), cex = cex.main, 
                    adj = 0.5, xpd = TRUE)
                } else {
                  text(x = lowerLeftCorner[1] + size[1] * 0.5, y = lowerLeftCorner[2] + 
                    size[2] * 1.05, labels = main[i], cex = cex.main, adj = 0.5)
                }
                
                if (pValue) {
                  
                  # check pValue slot
                  if (!is.null(testValue)) {
                    if (class(testValue) == "matrix") {
                      text(xat, apply(heightsOfRectangles, 2, max), signif(testValue[, 
                        i], 1), cex = cex.pValue, adj = c(0.5, 0.01))
                    } else {
                      text(xat, apply(heightsOfRectangles, 2, max), signif(testValue, 
                        1), cex = cex.pValue, adj = c(0.5, 0.01))
                    }
                    
                    
                  } else {
                    # use chisq by default
                    TFf <- fraction(x[i], strand = strand) < filter.pValue.fraction
                    ct <- signif(chisq.test(x[i], strand), 1)
                    ct[!TFf] <- NA
                    text(xat, apply(heightsOfRectangles, 2, max), ct, cex = cex.pValue, 
                      adj = c(0.5, 0.01))
                  }
                  
                }
                
                if (!is.null(OrgDb) | !is.null(TxDb)) {
					annotationBarplot(strand,i, lowerLeftCorner, annDfPlus, annDfMinus,
									  cex=cex.annotation, ypos=ypos.annotation,
									  interspace=annotation.interspace)
                }
                
                if (legend) {
                  
                  # any strand
                  o <- apply(height, 1, sum)
                  
                  # remove zeros count alleles
                  TFz <- rownames(height) == "X"
                  
                  o <- o[!TFz]
                  textOver <- rownames(height)[!TFz]
                  
                  # set another xLegendPos and yLegendPos for  '*'.
                  if ((strand == "*" ) & (!is.null(OrgDb) | 
                    !is.null(TxDb))) {
                    xLegendPos <- 0.5
                    yLegendPos <- 1.00
                  } else {
                    xLegendPos <- 0.94
                    yLegendPos <- 1.00
                  }
                  
				  legendBarplot(lowerLeftCorner, size, rownames=textOver, colnames=legend.colnames,
						   boxsize=legend.fill.size, boxspace=legend.interspace, fgCol, bgCol,
						   ylegendPos=1, xlegendPos=0.96, cex=cex.legend)

                  
                }
            }
            
        } else if (strand == "both") {
            for (i in 1:length(snps)) {
                
                if (verbose) {
                  cat("plotting snp", i, "\n")
                }
                
                snp <- snps[i]
                
                # getting revelant data forward strand
                tmpPlus <- t(alleleCounts(x, "+")[[snp]][samples, ])
                if ((ncol(tmpPlus) * nrow(tmpPlus)) == sum(tmpPlus == 0)) {
                  zeroRowsPlus <- TRUE
                } else {
                  zeroRowsPlus <- FALSE
                }
                
                tmpMinus <- t(alleleCounts(x, "-")[[snp]][samples, ])
                if ((ncol(tmpMinus) * nrow(tmpMinus)) == sum(tmpMinus == 0)) {
                  zeroRowsMinus <- TRUE
                } else {
                  zeroRowsMinus <- FALSE
                }
                
                # check most common alleles on + and - (should be same)
                allVar <- data.frame(A = c(NA, NA), T = c(NA, NA), G = c(NA, NA), 
                  C = c(NA, NA), del = c(NA, NA))
                o <- apply(tmpPlus, 1, sum)
                u <- apply(tmpMinus, 1, sum)
                
                # to get their total count
                allVar[1, match(names(o), colnames(allVar))] <- o
                allVar[2, match(names(u), colnames(allVar))] <- u
                
                # save to tri-allelic warner
                nAll <- names(sort(apply(allVar, 2, sum), decreasing = TRUE))
                tcAll <- sort(apply(allVar, 2, sum), decreasing = TRUE)
                
                # set colnames to X for non present alleles
                xNames <- colnames(allVar)[apply(allVar, 2, sum) == 0]
                colnames(allVar)[colnames(allVar) %in% xNames] <- "X"
                xNames <- rownames(tmpPlus)[apply(tmpPlus, 1, sum) == 0]
                rownames(tmpPlus)[rownames(tmpPlus) %in% xNames] <- "X"
                xNames <- rownames(tmpMinus)[apply(tmpMinus, 1, sum) == 0]
                rownames(tmpMinus)[rownames(tmpMinus) %in% xNames] <- "X"
                
                # alleles are
                n <- names(sort(apply(allVar, 2, sum), decreasing = TRUE))[1:2]
                # theor total count
                tc <- sort(apply(allVar, 2, sum), decreasing = TRUE)[1:2]
                
                # extracting only major and minor allele in correctly transposed format and
                # ordering colours accordingly
                over <- tmpPlus[match(n, rownames(tmpPlus)), ]
                under <- tmpMinus[match(n, rownames(tmpMinus)), ]
                # change all NAs to X names and zero values
                over[is.na(rownames(over)), ] <- 0
                rownames(over)[is.na(rownames(over))] <- "X"
                
                under[is.na(rownames(under)), ] <- 0
                rownames(under)[is.na(rownames(under))] <- "X"
                
                # Plus Warn if the two remaining alleles have too many counts
                if (!zeroRowsPlus) {
                  majorAndMinorFraction <- sum(allVar[1, ][c(n[1], n[2])])/sum(allVar[1, 
                    ])
                  if (verbose & majorAndMinorFraction < 0.9) {
                    message(paste("Snp", snp, "was possible tri-allelic, but only two most frequent alleles were plotted. Counts:\n", 
                      paste(paste(nAll, tcAll, sep = "="), collapse = ", ")))
                  }
                }
                # Minus calculating major and minor allele, warn if the two remaining alleles
                # have too many counts
                if (!zeroRowsMinus) {
                  majorAndMinorFraction <- sum(allVar[2, ][c(n[1], n[2])])/sum(allVar[2, 
                    ])
                  if (verbose & majorAndMinorFraction < 0.9) {
                    message(paste("Snp", snp, "was possible tri-allelic, but only two most frequent alleles were plotted. Counts:\n", 
                      paste(paste(nAll, tcAll, sep = "="), collapse = ", ")))
                  }
                }
                # set under as minus
                under <- -under
                
                # plotting function
                if (!add) {
                  # init frame
                  plot.default(NULL, xlim = c(0, 1), ylim = c(0, 1), ylab = "", xlab = "", 
                    xaxt = "n", yaxt = "n", , frame.plot = FALSE)
                }
                
                heightPlus <- over
                heightMinus <- under
                
                if (is.null(ylim)) {
                  if (zeroRowsPlus) {
                    ymax <- 5
                  } else {
                    ymax <- max(abs(range(heightPlus, na.rm = TRUE)))
                  }
                  if (zeroRowsMinus) {
                    ymin <- -5
                  } else {
                    ymin <- min(heightMinus)
                  }
                  ylimHere <- c(ymin, ymax)
                } else {
                  ylimHere <- ylim
                }
                
                # increase ylimHere before making pretty
                extraspaceOnTopAndBottom <- 0.15
                extraRealDistance <- (dist(ylimHere) * extraspaceOnTopAndBottom)
                ylimHere[1] <- ylimHere[1] - extraRealDistance
                ylimHere[2] <- ylimHere[2] + extraRealDistance
                
                # set lowest max and min values
                if (ylimHere[2] < 10) {
                  ylimHere[2] <- 10
                }
                if (abs(ylimHere[1]) < 10) {
                  ylimHere[1] <- -10
                }
                
                # make more even limits (make more pretty)
                ylabels <- pretty(ylimHere[1]:ylimHere[2])
                ylimHere[1] <- min(ylabels)
                ylimHere[2] <- max(ylabels)
                
                # breakpoint Plus Minus bars
                middleLine <- (1/(ylimHere[2] - ylimHere[1])) * abs(ylimHere[1])
                
                # make a xPointsAll framework
                xPointsAll <- seq(lowerLeftCorner[1], lowerLeftCorner[1] + size[1], 
                  length.out = length(heightPlus) * 1.5 + 2)[2:(length(heightPlus) * 
                  1.5 + 1)]
                
                widthsOfRectangles <- rep(size[1]/(length(heightPlus) * 1.5 + 1), 
                  length(heightPlus)/2)
                # normalize toPlot-values to a range of c(0,size[2])
                heightsOfRectanglesPlus <- round(((heightPlus - ylimHere[1])/(ylimHere[2] - 
                  ylimHere[1]) * size[2]) - middleLine, 5)
                heightsOfRectanglesMinus <- round(((-heightMinus - ylimHere[1])/(ylimHere[2] - 
                  ylimHere[1]) * size[2]) - middleLine, 5)
                
                # take every third element (start at 2)
                xPoints1 <- xPointsAll[seq(2, length(heightPlus) * 1.5, by = 3)]
                yPointsPlus <- lowerLeftCorner[2] + middleLine + heightsOfRectanglesPlus[1, 
                  ]/2
                
                barMinusEdge <- matrix(rep(NA, length(heightsOfRectanglesMinus)), 
                  nrow = 2)  #is used later to plot p-value
                
                # yPointMinus is a little tricky and needs a more complicated equation
                yPointsMinus1 <- lowerLeftCorner[2] + heightsOfRectanglesMinus[1, 
                  ]/2
                yPointsMinus2 <- 1 - yPointsMinus1
                yPointsMinus3 <- yPointsMinus2 - (1 - middleLine)
                
                barMinusEdge[1, ] <- (1 - (lowerLeftCorner[2] + heightsOfRectanglesMinus[1, 
                  ])) - (1 - middleLine)
                
                yPointsMinus <- yPointsMinus3
                
                alleleColors1 <- fgCol
                
                symbols(x = xPoints1, y = yPointsPlus, bg = alleleColors1, rectangles = matrix(ncol = 2, 
                  c(widthsOfRectangles, heightsOfRectanglesPlus[1, ])), add = TRUE, 
                  inches = FALSE)
                symbols(x = xPoints1, y = yPointsMinus, bg = alleleColors1, rectangles = matrix(ncol = 2, 
                  c(widthsOfRectangles, heightsOfRectanglesMinus[1, ])), add = TRUE, 
                  inches = FALSE)
                
                # take every third element (start at 3)
                xPoints2 <- xPointsAll[seq(3, length(heightPlus) * 1.5, by = 3)]
                yPointsPlus <- lowerLeftCorner[2] + middleLine + heightsOfRectanglesPlus[2, 
                  ]/2
                
                # yPointMinus is a little tricky and needs a more complicated equation
                yPointsMinus1 <- lowerLeftCorner[2] + heightsOfRectanglesMinus[2, 
                  ]/2
                yPointsMinus2 <- 1 - yPointsMinus1
                yPointsMinus3 <- yPointsMinus2 - (1 - middleLine)
                
                barMinusEdge[2, ] <- (1 - (lowerLeftCorner[2] + heightsOfRectanglesMinus[2, 
                  ])) - (1 - middleLine)
                
                yPointsMinus <- yPointsMinus3
                
                # alleleColors2 <- rep('#7F5959',9)
                alleleColors2 <- bgCol
                
                symbols(x = xPoints2, y = yPointsPlus, bg = alleleColors2, rectangles = matrix(ncol = 2, 
                  c(widthsOfRectangles, heightsOfRectanglesPlus[2, ])), add = TRUE, 
                  inches = FALSE)
                symbols(x = xPoints2, y = yPointsMinus, bg = alleleColors2, rectangles = matrix(ncol = 2, 
                  c(widthsOfRectangles, heightsOfRectanglesMinus[2, ])), add = TRUE, 
                  inches = FALSE)
                
                # y-axis
                yat <- (ylabels + abs(min(ylabels)))/(sum(abs(min(ylabels)) + max(ylabels)))
                graphParamList[["yat"]] <- yat
                if (yaxis) {
                  # vertical line
                  lines(c(lowerLeftCorner[1], lowerLeftCorner[1]), c(lowerLeftCorner[2], 
                    lowerLeftCorner[2] + 1))
                  # horizontal lines
                  invisible(mapply(FUN = function(x, y) {
                    lines(c(x, x - size[1] * 0.01), c(y, y))
                  }, rep(lowerLeftCorner[1], length(yat)), (lowerLeftCorner[2] + 
                    yat)))
                }
                if (ylab) {
                  if (las.ylab == 2) {
                    text(lowerLeftCorner[1] - size[1] * 0.02, yat, labels = abs(ylabels), 
                      xpd = TRUE, srt = 90, adj = c(1, 0.5), cex = cex.ylab)
                  }
                  if (las.ylab == 1) {
                    text(lowerLeftCorner[1] - size[1] * 0.02, yat, labels = abs(ylabels), 
                      xpd = TRUE, srt = 0, adj = c(1, 0.5), cex = cex.ylab)
                  }
                }
                
                # x-axis labels
                xat <- xPoints1 + (xPoints2 - xPoints1)/2
                graphParamList[["xat"]] <- xat
                if (xaxis) {
                  # vertical line
                  lines(c(lowerLeftCorner[1], lowerLeftCorner[1] + 1), c(lowerLeftCorner[2], 
                    lowerLeftCorner[2]))
                  # horizontal lines
                  invisible(mapply(FUN = function(x, y) {
                    lines(c(x, x), c(y, y - size[1] * 0.01))
                  }, (lowerLeftCorner[2] + xat), rep(lowerLeftCorner[1], length(xat))))
                }
                if (xlab) {
                  if (las.xlab == 2) {
                    text(xat, lowerLeftCorner[2] - size[2] * 0.02, labels = samples, 
                      xpd = TRUE, srt = 90, adj = c(1, 0.5), cex = cex.xlab)
                  }
                  if (las.xlab == 1) {
                    text(xat, lowerLeftCorner[2] - size[2] * 0.02, labels = samples, 
                      xpd = TRUE, srt = 0, adj = c(1, 0.5), cex = cex.xlab)
                  }
                }
                if (is.null(main)) {
                  text(x = lowerLeftCorner[1] + size[1] * 0.5, y = lowerLeftCorner[2] + 
                    size[2] * 1.05, labels = paste(snp), cex = cex.main, 
                    adj = 0.5, xpd = TRUE)
                } else {
                  text(x = lowerLeftCorner[1] + size[1] * 0.5, y = lowerLeftCorner[2] + 
                    size[2] * 1.05, labels = main[i], cex = cex.main, adj = 0.5)
                }
                
                
                
                if (pValue) {
                  
                  # check pValue slot
                  if (!is.null(testValue)) {
                    if (class(testValue) == "matrix") {
                      text(xat, apply(middleLine + heightsOfRectanglesPlus, 2, max), 
                        signif(testValue[, i], 1), cex = cex.pValue, adj = c(0.4, 
                          -0.3))
                    } else {
                      text(xat, apply(middleLine + heightsOfRectanglesPlus, 2, max), 
                        signif(testValue, 1), cex = cex.pValue, adj = c(0.4, -0.3))
                    }
                  } else {
                    # use chisq by default
                    TFf <- fraction(x[i], strand = "+") < filter.pValue.fraction
                    ct <- signif(chisq.test(x[i], "+"), 1)
                    ct[!TFf] <- NA
                    text(xat, apply(middleLine + heightsOfRectanglesPlus, 2, max), 
                      ct, cex = cex.pValue, adj = c(0.4, -0.3))
                    
                  }
                  if (!is.null(testValue2)) {
                    if (class(testValue) == "matrix") {
                      text(xat, apply(barMinusEdge, 2, min), signif(testValue2[, 
                        i], 1), cex = cex.pValue, adj = c(0.4, 1, 1))
                    } else {
                      text(xat, apply(barMinusEdge, 2, min), signif(testValue2, 1), 
                        cex = cex.pValue, adj = c(0.4, 1, 1))
                    }
                  } else {
                    # use chisq by default
                    TFf <- fraction(x[i], strand = "-") < filter.pValue.fraction
                    ct <- signif(chisq.test(x[i], "-"), 1)
                    ct[!TFf] <- NA
                    text(xat, apply(barMinusEdge, 2, min), ct, cex = cex.pValue, 
                      adj = c(0.4, 1, 1))
                    
                  }
                  
                }
                
                if (!is.null(OrgDb) | !is.null(TxDb)) {
                  
					annotationBarplot(strand,i, lowerLeftCorner, annDfPlus, annDfMinus,
									  cex=cex.annotation, ypos=ypos.annotation,
									  interspace=annotation.interspace)
#                  # check which columns to extract
#                  dfType <- c("symbol", "exon_id", "tx_id", "cds_id")
#                  TFann <- (dfType %in% colnames(annDfPlus)) | (dfType %in% colnames(annDfMinus))
#                  
#                  # plus strand (over)
#                  df <- annDfPlus[, TFann, drop = FALSE]
#                  yPosTmp <- lowerLeftCorner[2] + 1.01
#                  for (j in 1:length(colnames(df))) {
#                    text <- paste(colnames(df)[j], df[i, j], sep = ": ")
#                    text(x = lowerLeftCorner[1] + 0.01, y = yPosTmp, labels = text, 
#                      cex = 0.7, adj = 0)
#                    yPosTmp <- yPosTmp - 0.02
#                  }
#                  
#                  # minus strand (under)
#                  df <- annDfMinus[, TFann, drop = FALSE]
#                  
#                  yPosTmp <- lowerLeftCorner[2] - 0.01
#                  for (j in 1:length(colnames(df))) {
#                    text <- paste(colnames(df)[j], df[i, j], sep = ": ")
#                    text(x = lowerLeftCorner[1] + 0.01, y = yPosTmp, labels = text, 
#                      cex = 0.7, adj = 0)
#                    yPosTmp <- yPosTmp + 0.02
#                  }
                }
                
                if (legend) {
                  # plus strand
                  TFz <- rownames(over) == "X"

                  if (!sum(!TFz) == 0) {
                    textOver <- paste(rownames(over)[!TFz], apply(over, 1, sum)[!TFz])
					#remove the total count for that nucleotide
					textOver <- unlist(lapply(strsplit(textOver," "),function(x){x[1]}))

					  legendBarplot(lowerLeftCorner, size, rownames=textOver, colnames=legend.colnames,
						   boxsize=legend.fill.size, boxspace=legend.interspace, fgCol, bgCol,
						   ylegendPos=1, xlegendPos=0.96, cex=cex.legend)
                    }
                  }
                 # # minus strand
                 # TFz <- rownames(under) == "X"
                 # if (!sum(!TFz) == 0) {
                 #   textUnder <- paste(rownames(under)[!TFz], -apply(under, 1, sum)[!TFz])
                 #   
				 #   #remove the total count for that nucleotide
				 #   textUnder <- unlist(lapply(strsplit(textUnder," "),function(x){x[1]}))

				 #   #user option for changing the legend.size
				 #   if(!is.null(legend.fill.size)){
				 #   	size <- legend.fill.size
				 #   }else{
				 #   	size <- size
				 #   }

                 #   if (sum(!TFz) == 2) {
                 #     symbols(x = lowerLeftCorner[1] + size[1] * seq(0.96, (0.96 - 
                 #       (0.02 * (length(unique(fgCol)) - 1))), length = length(unique(fgCol))), 
                 #       y = (lowerLeftCorner[2] + size[2] * 0.03) * rep(1, length(unique(fgCol))), 
                 #       bg = unique(fgCol), squares = rep(c(size[1] * 0.012), length(unique(fgCol))), 
                 #       add = TRUE, inches = FALSE)
                 #     symbols(x = lowerLeftCorner[1] + size[1] * seq(0.96, (0.96 - 
                 #       (0.02 * (length(unique(bgCol)) - 1))), length = length(unique(bgCol))), 
                 #       y = (lowerLeftCorner[2] + size[2] * 0.01) * rep(1, length(unique(bgCol))), 
                 #       bg = unique(bgCol), squares = rep(c(size[1] * 0.012), length(unique(bgCol))), 
                 #       add = TRUE, inches = FALSE)
                 #     text(x = c(lowerLeftCorner[1] + (size[1] * c(0.98, 0.98))), 
                 #       y = lowerLeftCorner[2] + size[2] * c(0.03, 0.01), textUnder, 
                 #       srt = 0, cex = cex.legend, adj = c(0, 0.5), xpd = TRUE)
                 #   } else if (sum(!TFz[1]) == 1) {
                 #     symbols(x = lowerLeftCorner[1] + size[1] * seq(0.96, (0.96 - 
                 #       (0.02 * (length(unique(fgCol)) - 1))), length = length(unique(fgCol))), 
                 #       y = (lowerLeftCorner[2] + size[2] * 0.01) * rep(1, length(unique(fgCol))), 
                 #       bg = unique(fgCol), squares = rep(c(size[1] * 0.012), length(unique(fgCol))), 
                 #       add = TRUE, inches = FALSE)
                 #     text(x = c(lowerLeftCorner[1] + (size[1] * c(0.98))), y = lowerLeftCorner[2] + 
                 #       size[2] * c(0.01), textUnder, srt = 0, cex = cex.legend, 
                 #       adj = c(0, 0.5), xpd = TRUE)
                 #     
                 #   } else if (sum(!TFz[2]) == 1) {
                 #     symbols(x = lowerLeftCorner[1] + size[1] * seq(0.96, (0.96 - 
                 #       (0.02 * (length(unique(bgCol)) - 1))), length = length(unique(bgCol))), 
                 #       y = (lowerLeftCorner[2] + size[2] * 0.01) * rep(1, length(unique(bgCol))), 
                 #       bg = unique(bgCol), squares = rep(c(size[1] * 0.012), length(unique(bgCol))), 
                 #       add = TRUE, inches = FALSE)
                 #     text(x = c(lowerLeftCorner[1] + (size[1] * c(0.98, 0.98))), 
                 #       y = lowerLeftCorner[2] + size[2] * c(0.01), textUnder, srt = 0, 
                 #       cex = cex.legend, adj = c(0, 0.5), xpd = TRUE)
                 #   }
                 #   
                 # }
                #}
            }
        }
    #} else if (type == "fraction") {
    #    # a plot type that shows the allelic imbalance, irrespective of absolute read
    #    # counts at each Snp.  this produce plots that are more easily overviewed than
    #    # those from 'plot_reads_by_allele', but of course omits some detail information
    #    # from these plots. Intended to use for a medium number of snps and samples
    #    
    #    if (strand == "both") {
    #        stop("strand must be '+', '-' or '*' for type='fraction' ")
    #    }
    #    
    #    for (i in 1:length(snps)) {
    #        # getting revelant data
    #        snp <- snps[i]
    #        fgColHere <- fgCol
    #        tmp <- alleleCounts(x, strand)[[snp]][samples, , drop = FALSE]
    #        
	#		#df <- fractionPlotDf(e$x, identifier, strand=e$strand, top.allele.criteria=e$top.allele.criteria)
    #        # check for zeroRows
    #        if ((ncol(tmp) * nrow(tmp)) == sum(tmp == 0)) {
    #            zeroRows <- TRUE
    #        } else {
    #            zeroRows <- FALSE
    #        }
    #        
    #        # calculating major and minor allele, warn if the two remaining alleles have too
    #        # many counts
    #        if (nrow(tmp) > 1) {
    #            countsByAllele <- apply(tmp, 2, sum, na.rm = TRUE)
    #        } else {
    #            countsByAllele <- tmp[1, ]
    #        }
    #        
    #        if (!zeroRows) {
    #            
    #            majorAllele <- colnames(tmp)[order(countsByAllele, decreasing = TRUE)][1]
    #            minorAllele <- colnames(tmp)[order(countsByAllele, decreasing = TRUE)][2]
    #            majorAndMinorFraction <- sum(countsByAllele[c(majorAllele, minorAllele)])/sum(countsByAllele)
    #            if (verbose & majorAndMinorFraction < 0.9) {
    #              message(paste("Snp", snp, "was possible tri-allelic, but only two most frequent alleles were plotted. Counts:"))
    #              message(paste(paste(names(countsByAllele), countsByAllele, sep = "="), 
    #                collapse = ", "))
    #            }
    #            
    #            fraction <- tmp[, majorAllele]/(tmp[, minorAllele] + tmp[, majorAllele])
    #        } else {
    #            fraction <- tmp[, 1]/(tmp[, 1] + tmp[, 2])
    #        }
    #        
    #        # setting no-count samples to colour grey90 and fraction 0
    #        fgColHere[is.nan(fraction)] <- "grey90"
    #        fraction[is.nan(fraction)] <- 0
    #        
    #        # if add=FALSE then make a new device
    #        if (!add) {
    #            plot.default(NULL, xlim = c(0, 1), ylim = c(0, 1), ylab = "", xlab = "", 
    #              xaxt = "n", yaxt = "n")
    #        }
    #        
    #        # ylim is always (0,1) and usergiven ylims are reset to this with a warning (in
    #        # start checkups).
    #        ylim <- c(0, 1)
    #        
    #        # find the centers of rectangles
    #        xPoints <- seq(lowerLeftCorner[1], lowerLeftCorner[1] + size[1], length.out = length(fraction) + 
    #            2)[2:(length(fraction) + 1)]
    #        widthsOfRectangles <- rep(size[1]/(length(fraction) + 1), length(fraction))
    #        
    #        # normalize fraction-values to a range of c(0,size[2])
    #        heightsOfRectangles <- ((fraction - ylim[1])/(ylim[2] - ylim[1])) * size[2]
    #        # find the Y-coordinate center points for rectangles
    #        yPoints <- lowerLeftCorner[2] + heightsOfRectangles/2
    #        
    #        # if bgCol is given, start by painting background
    #        if (!is.null(bgCol)) {
    #            yPointsBg <- rep(lowerLeftCorner[2] + size[2]/2, length(xPoints))
    #            rectanglesBg <- matrix(ncol = 2, c(widthsOfRectangles, rep(size[2], 
    #              length(widthsOfRectangles))))
    #            symbols(x = xPoints, y = yPointsBg, bg = bgCol, rectangles = rectanglesBg, 
    #              add = TRUE, inches = FALSE)
    #        }
    #        
    #        symbols(x = xPoints, y = yPoints, bg = fgColHere, rectangles = matrix(ncol = 2, 
    #            c(widthsOfRectangles, heightsOfRectangles)), add = TRUE, inches = FALSE)
    #        
    #        # add a frame to the plot
    #        if (add.frame) {
    #            symbols(x = lowerLeftCorner[1] + size[1]/2, y = lowerLeftCorner[2] + 
    #              size[2]/2, rectangles = matrix(ncol = 2, c(size[1], size[2])), 
    #              add = TRUE, inches = FALSE)
    #        }
    #        
    #        # add a tick markers (easier for fraction-plots because they are restricted to
    #        # 0-100%
    #        marks <- c(0.25, 0.5, 0.75)
    #        names(marks) <- c("25%", "50%", "75%")
    #        for (markName in names(marks)) {
    #            mark <- marks[markName]
    #            yBase <- lowerLeftCorner[2] + ((mark - ylim[1])/(ylim[2] - ylim[1])) * 
    #              size[2]
    #            xBase <- xPoints[1] - widthsOfRectangles[1]/2
    #            lines(x = c(xBase - size[1] * 0.01, xBase + size[1] * 0.01), y = c(yBase, 
    #              yBase))
    #            text(x = xBase - size[1] * 0.015, y = yBase, label = markName, adj = 1, 
    #              cex = 0.7)
    #        }
    #        
    #        # writing the name of each Snp on top of the plot(within frame)
    #        text(x = 0.5, y = 1.02, label = snp, xpd = TRUE)
    #        
    #        # add horizontal line at 0.5
    #        if (!is.null(addHorizontalLine)) {
    #            lines(x = c(lowerLeftCorner[1], lowerLeftCorner[1] + size[1]), y = c(lowerLeftCorner[2] + 
    #              size[2] * addHorizontalLine, lowerLeftCorner[2] + size[2] * addHorizontalLine))
    #            
    #            mapBiasHere <- mapBias(x)[[snp]]
    #            if (!all(mapBiasHere %in% c(0.5, 0)) & !zeroRows) {
    #              # Only activate if anything is different from 0.5 or 0
    #              
    #              for (j in 1:length(fraction)) {
    #                biasfraction <- mapBiasHere[j, majorAllele]/(mapBiasHere[j, majorAllele] + 
    #                  mapBiasHere[j, minorAllele])
    #                yPoint <- lowerLeftCorner[2] + ((biasfraction - ylim[1])/(ylim[2] - 
    #                  ylim[1])) * size[2]
    #                lines(x = c(xPoints[j] - widthsOfRectangles[j]/2, xPoints[j] + 
    #                  widthsOfRectangles[j]/2), y = c(yPoint, yPoint), lty = 2)
    #              }
    #            }
    #        }
    #    }
    #}
    } else if (type == "fraction") {
        # a plot type that shows the allelic imbalance, irrespective of absolute read
        # counts at each Snp.  this produce plots that are more easily overviewed than
        # those from 'plot_reads_by_allele', but of course omits some detail information
        # from these plots. Intended to use for a medium number of snps and samples
        
        if (strand == "both") {
            stop("strand must be '+', '-' or '*' for type='fraction' ")
        }
        
        for (i in 1:length(snps)) {
            # getting revelant data
            snp <- snps[i]
            fgColHere <- fgCol
            #tmp <- alleleCounts(x, strand)[[snp]][samples, , drop = FALSE]
            
			df <- fractionPlotDf(x, snp, strand=strand, top.allele.criteria=top.allele.criteria)
            # check for zeroRows
            if (sum(df$na=="no")==0) {
                zeroRows <- TRUE
            } else {
                zeroRows <- FALSE
            }
            
            # setting no-count samples to colour grey90 and fraction 0
            #fgColHere[df$na[seq(1,nrow(df),by=2)+1]=="yes"] <- "grey90"
            bgCol[df$na[seq(1,nrow(df),by=2)+1]=="yes"] <- "grey90"
            #fraction[is.nan(fraction)] <- 0
            
            # if add=FALSE then make a new device
            if (!add) {
                plot.default(NULL, xlim = c(0, 1), ylim = c(0, 1), ylab = "", xlab = "", 
                  xaxt = "n", yaxt = "n")
            }
            
            # ylim is always (0,1) and usergiven ylims are reset to this with a warning (in
            # start checkups).
            ylim <- c(0, 1)
            
            # find the centers of rectangles
            xPoints <- seq(lowerLeftCorner[1], lowerLeftCorner[1] + size[1], length.out = nrow(df)/2 + 
                2)[2:(nrow(df)/2 + 1)]
            widthsOfRectangles <- rep(size[1]/(nrow(df)/2 + 1), nrow(df)/2)
            
            # normalize fraction-values to a range of c(0,size[2])
            heightsOfRectangles <- ((df$values[seq(1,nrow(df),by=2)] - ylim[1])/(ylim[2] - ylim[1])) * size[2]
            # find the Y-coordinate center points for rectangles
            yPoints <- lowerLeftCorner[2] + heightsOfRectangles/2
            
            # if bgCol is given, start by painting background
            if (!is.null(bgCol)) {
                yPointsBg <- rep(lowerLeftCorner[2] + size[2]/2, length(xPoints))
                rectanglesBg <- matrix(ncol = 2, c(widthsOfRectangles, rep(size[2], 
                  length(widthsOfRectangles))))
                symbols(x = xPoints, y = yPointsBg, bg = bgCol, rectangles = rectanglesBg, 
                  add = TRUE, inches = FALSE)
            }
            
            symbols(x = xPoints, y = yPoints, bg = fgColHere, rectangles = matrix(ncol = 2, 
                c(widthsOfRectangles, heightsOfRectangles)), add = TRUE, inches = FALSE)
            
            # add a frame to the plot
            if (add.frame) {
                symbols(x = lowerLeftCorner[1] + size[1]/2, y = lowerLeftCorner[2] + 
                  size[2]/2, rectangles = matrix(ncol = 2, c(size[1], size[2])), 
                  add = TRUE, inches = FALSE)
            }
            
            # add a tick markers (easier for fraction-plots because they are restricted to
            # 0-100%
            marks <- c(0.25, 0.5, 0.75)
            names(marks) <- c("25%", "50%", "75%")
            for (markName in names(marks)) {
                mark <- marks[markName]
                yBase <- lowerLeftCorner[2] + ((mark - ylim[1])/(ylim[2] - ylim[1])) * 
                  size[2]
                xBase <- xPoints[1] - widthsOfRectangles[1]/2
                lines(x = c(xBase - size[1] * 0.01, xBase + size[1] * 0.01), y = c(yBase, 
                  yBase))
                text(x = xBase - size[1] * 0.015, y = yBase, label = markName, adj = 1, 
                  cex = 0.7)
            }
            
            # writing the name of each Snp on top of the plot(within frame)
            text(x = 0.5, y = 1.02, label = snp, xpd = TRUE, cex = cex.main)
            
            # add horizontal line at 0.5
            if (!is.null(addHorizontalLine)) {
                lines(x = c(lowerLeftCorner[1], lowerLeftCorner[1] + size[1]), y = c(lowerLeftCorner[2] + 
                  size[2] * addHorizontalLine, lowerLeftCorner[2] + size[2] * addHorizontalLine))
                
                mapBiasHere <- mapBias(x)[[snp]]
                #if (!all(mapBiasHere %in% c(0.5, 0)) & !zeroRows) {
                #  # Only activate if anything is different from 0.5 or 0
                #  
                #  for (j in 1:length(fraction)) { #fraction doesnt exist (check df$values instead)
                #    biasfraction <- mapBiasHere[j, majorAllele]/(mapBiasHere[j, majorAllele] + 
                #      mapBiasHere[j, minorAllele])
                #    yPoint <- lowerLeftCorner[2] + ((biasfraction - ylim[1])/(ylim[2] - 
                #      ylim[1])) * size[2]
                #    lines(x = c(xPoints[j] - widthsOfRectangles[j]/2, xPoints[j] + 
                #      widthsOfRectangles[j]/2), y = c(yPoint, yPoint), lty = 2)
                #  }
                #}
            }
        }
    }
    # return graphical paramlist (does not work for fraction type yet)
    invisible(graphParamList)
})

#' lbarplot ASEset objects
#' 
#' Generates lbarplots for ASEset objects. Two levels of plotting detail are
#' provided: a detailed lbarplot of read counts by allele useful for fewer
#' samples and SNPs, and a less detailed lbarplot of the fraction of imbalance,
#' useful for more samples and SNPs.
#' 
#' This function serves the same purpose as the normal barplot, but with
#' trellis graphics using lattice, to be able to integrate well with Gviz track
#' functionality.
#' 
#' @name ASEset-lbarplot
#' @aliases ASEset-lbarplot lbarplot lbarplot,ASEset-method
#' @docType methods
#' @param x An \code{ASEset} object
#' @param type 'count' or 'fraction'
#' @param strand four options, '+', '-', 'both' or '*'
#' @param verbose Makes function more talkative
#' @param ... for simpler generics when extending function
#' @author Jesper R. Gadin
#' @seealso \itemize{ \item The \code{\link{ASEset}} class which the lbarplot
#' function can be called up on.  \item The \code{\link{barplot}} non trellis
#' barplot.  }
#' @keywords lbarplot
#' @examples
#' 
#' data(ASEset)
#' lbarplot(ASEset[1])
#' 
#' @exportMethod lbarplot

setGeneric("lbarplot", function(x, type = "count", strand = "*", 
    verbose = FALSE, ...) {
    standardGeneric("lbarplot")
})

setMethod("lbarplot", signature(x = "ASEset"), function(x, type = "count", strand = "*", 
    usePhase=FALSE, verbose = FALSE, ...) {

    if (length(list(...)) == 0) {
        e <- new.env(hash = TRUE)
    } else {
        e <- list2env(list(...))
    }

	#print(ls(envir=e))

    if (!exists("mainvec", envir = e, inherits = FALSE)) {
		e$mainvec <- rep("",nrow(x))
	}
    if (!exists("cex.mainvec", envir = e, inherits = FALSE)) {
		e$cex.mainvec <- 1
	}
    if (!exists("ylab", envir = e, inherits = FALSE)) {
        e$ylab <- ""
    }
    if (!exists("xlab", envir = e, inherits = FALSE)) {
        e$xlab <- ""
    }
    if (!exists("middleLine", envir = e, inherits = FALSE)) {
        e$middleLine <- TRUE
    }
    if (!exists("deAnnoPlot", envir = e, inherits = FALSE)) {
        e$deAnnoPlot <- FALSE
    }
    if (!exists("top.allele.criteria", envir = e, inherits = FALSE)) {
        e$top.allele.criteria <- "maxcount"
    }

	for (i in 1:nrow(x)) {
		name <- rownames(x)[i]
		if (type == "fraction") {
			b <- barplotLatticeFraction(
				identifier = name, 
				x=x, 
				astrand=strand, 
				ids=rownames(x),
				ylab=e$ylab,
				xlab=e$xlab,
				mainvec=e$mainvec,
				middleLine=e$middleLine,
				deAnnoPlot=e$deAnnoPlot,
				cex.mainvec=e$cex.mainvec,
				top.allele.criteria=e$top.allele.criteria)
		} else if (type == "count") {
			b <- barplotLatticeCounts(
				identifier = name, 
				x=x, 
				astrand=strand, 
				ids=rownames(x),
				ylab=e$ylab,
				xlab=e$xlab,
				mainvec=e$mainvec)
		} else {
			stop("type has to be fraction or count")
		}
	}
    b
}) 

