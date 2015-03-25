#'@include ASEset-class.R
NULL





#' Get Gene Area
#' 
#' Given a character vector with genesymbols and an OrgDb object, this function
#' returns a GRanges giving the coordinates of the genes.
#' 
#' This function is a convenience function that can be used to determine which
#' genomic coordinates to specify to e.g. \code{impBamGAL} when retrieving
#' reads.
#' 
#' The function cannot handle genes that do not exist in the annotation. To
#' remove these please set the na.rm=TRUE.
#' 
#' @param genesymbols A character vector that contains genesymbols of genes
#' from which we wish to retrieve the coordinates
#' @param OrgDb An \code{OrgDb} object containing gene annotation
#' @param leftFlank A \code{integer} specifying number of additional
#' nucleotides before the genes
#' @param rightFlank A \code{integer} specifying number of additional
#' nucleotides after the genes
#' @param na.rm A \code{boolean} removing genes that returned NA from the
#' annotation
#' @param verbose Setting \code{verbose=TRUE} makes function more talkative
#' @return \code{getAreaFromGeneNames} returns a GRanges object with genomic
#' coordinates around the specified genes
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords genes locations
#' @examples
#' 
#' #load example data
#' data(ASEset)
#' 
#' #get counts at the three positions specified in GRvariants
#' library(org.Hs.eg.db )
#' searchArea<-getAreaFromGeneNames(c('PAX8','TLR7'), org.Hs.eg.db)
#' 
#' @export getAreaFromGeneNames
getAreaFromGeneNames <- function(genesymbols, OrgDb, leftFlank = 0, rightFlank = 0, 
    na.rm = FALSE, verbose = TRUE) {
    
    # start up sets
    if (!class(genesymbols) %in% c("character")) 
        stop(paste("genesymbols must be of class character, not:", class(genesymbols)))
    
    if (!class(OrgDb) %in% c("OrgDb")) 
        stop(paste("OrgDb must be of class OrgDb, not:", class(OrgDb)))
    
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
    
    # retrieving data
    colsFilter <- c("CHR", "CHRLOC", "CHRLOCEND", "SYMBOL")
    s <- suppressWarnings(select(OrgDb, keys = genesymbols, cols = colsFilter, keytype = "SYMBOL"))
    
    missing <- genesymbols[!genesymbols %in% s[, "SYMBOL"]]
    if (verbose & length(missing) > 0) {
        cat(paste("Did not find information on these", length(missing), "genes:", 
            paste(missing, collapse = ", ")), "\n")
    } else {
        if (verbose) 
            cat("Found all requested genes in annotation", "\n")
    }
    
    # remove NAs
    if (na.rm == TRUE) {
        s <- s[!(is.na(s[, "CHR"]) | is.na(s[, "CHRLOC"]) | is.na(s[, "CHRLOCEND"])), 
            ]
    } else {
        warningGenes <- s[is.na(s[, "CHR"]) | is.na(s[, "CHRLOC"]) | is.na(s[, "CHRLOCEND"]), 
            ]
        if (!nrow(warningGenes) == 0) {
            warning(paste(warningGenes[, "SYMBOL"], "had NAs", "\n"), "you better remove these genes from your 'genesymbols'")
        }
    }
    
    strand <- rep("+", nrow(s))
    TFstrand <- s[, "CHRLOC"] < 0
    strand[TFstrand] <- "-"
    
    searchArea <- GRanges(seqnames = paste("chr", s[, "CHR"], sep = ""), ranges = IRanges(abs(s[, 
        "CHRLOC"]) - leftFlank, abs(s[, "CHRLOCEND"]) + rightFlank), strand = strand, 
        symbol = s[["SYMBOL"]])
    
    searchArea <- reduce(searchArea, with.revmap = TRUE)
    l <- lapply(searchArea$revmap, function(x) {
        paste(unique(s[x, "SYMBOL"]), collapse = ",")
    })
    mcols(searchArea)[, "symbol"] <- unlist(l)
    # remove the mapping column
    mcols(searchArea)[["revmap"]] <- NULL
    
    return(searchArea)
}


#' Get rsIDs from locations of SNP
#' 
#' Given a GRanges object of SNPs and a SNPlocs annotation, this function
#' attempts to replace the names of the GRanges object entries with rs-IDs.
#' 
#' This function is used to try to identify the rs-IDs of SNPs in a GRanges
#' object.
#' 
#' @param GR A \code{GRanges} that contains positions of SNPs to look up
#' @param SNPloc A \code{SNPlocs object} containing information on SNP
#' locations (e.g. SNPlocs.Hsapiens.dbSNP.xxxxxxxx)
#' @param return.vector Setting \code{return.vector=TRUE} returns vector with
#' rsIds
#' @param verbose Setting \code{verbose=TRUE} makes function more talkative
#' @return \code{getSnpIdFromLocation} returns the same GRanges object it was
#' given with, but with updated with rs.id information.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords SNP rs-id
#' @examples
#' 
#' #load example data
#' data(ASEset)
#' 
#' #get counts at the three positions specified in GRvariants
#' if(require(SNPlocs.Hsapiens.dbSNP.20120608)){
#' updatedGRanges<-getSnpIdFromLocation(rowRanges(ASEset), SNPlocs.Hsapiens.dbSNP.20120608)
#' rowRanges(ASEset)<-updatedGRanges
#'}
#' 
#' @export getSnpIdFromLocation
getSnpIdFromLocation <- function(GR, SNPloc, return.vector = FALSE, verbose = TRUE) {
   
	#genome <- genome(GR)

    if (class(GR) != "GRanges") 
        stop(paste("GR must of class GRanges, not", class(GR)))
    if (class(SNPloc) != "SNPlocs") 
        stop(paste("SNPlocs must of class SNPlocs, not", class(SNPloc)))
    if (!exists("getSNPlocs")) 
        stop("There must exist a function called getSNPlocs, available from the SNPlocs.Hsapiens.dbSNP.xxxxxxx package. Try to load such package.")
    
    # add chr to seqnames if not present
    if (length(grep("^chr", seqnames(GR))) != length(GR)) {
        if (length(grep("^chr", seqnames(GR))) != 0) 
            stop("seqnames must all begin with 'chr'. In the GR it seemed that some, but not all, seqnames began with chr. Please correct manually")
        seqlevels(GR) <- paste("chr", seqlevels(GR), sep = "")
        # seqnames(GR)<-seqnames
    }
    
    # changing chr to ch to adapt to SNPloc
    if (length(grep("^chr", seqnames(GR))) == length(GR)) {
        # seqnames<-seqnames(GR)
        seqlevels(GR) <- sub("^chr", "ch", seqlevels(GR))
        # seqlevels(GR)<-as.character(unique(seqnames)) seqlevels(GR) <- levels(seqnames)
        # seqnames(GR)<-seqnames
    }
    
    
    SNPlocThisChr <- getSNPlocs(seqlevels(GR), as.GRanges = TRUE, caching = FALSE)
    
    seqlevels(SNPlocThisChr, force = TRUE) <- seqlevels(GR)
    # seqlengths(GR) <- seqlengths(SNPlocThisChr)
	genome(SNPlocThisChr) <- genome(GR)
    
    overlaps <- findOverlaps(GR, SNPlocThisChr)
    
    if (verbose) 
        cat(paste("Replacing position-based SNP name with rs-ID for", length(overlaps), 
            "SNP(s)"), "\n")
    
    # replace name in GR for(i in 1:length(overlaps)){
    # snp<-paste('rs',mcols(SNPlocThisChr[subjectHits(overlaps[i])])[,'RefSNP_id'],sep='')
    # names(GR)[queryHits(overlaps[i])] <-snp }
    
    snp <- paste("rs", mcols(SNPlocThisChr[subjectHits(overlaps)])[, "RefSNP_id"], 
        sep = "")
    names(GR)[queryHits(overlaps)] <- snp
    
    # change back to chr from ch
    seqlevels(GR) <- sub("^ch", "chr", seqlevels(GR))


    if (return.vector) {
        names(GR)
    } else {
        return(GR)
    }
    
}

#' lattice barplot inner functions for ASEset objects
#' 
#' Generates lattice barplots for ASEset objects. Two levels of plotting detail
#' are provided: a detailed barplot of read counts by allele useful for fewer
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
#' @name barplot-lattice-support
#' @rdname barplot-lattice-support
#' @aliases barplot-lattice-support barplotLatticeCounts barplotLatticeFraction
#' @param identifier, the single snp name to plot
#' @param ... used to pass on variables
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{ASEset}} class which the barplot
#' function can be called up on.  }
#' @keywords barplot
#' @examples
#' 
#' a <- ASEset
#' name <- rownames(a)[1]
#' 
#' barplotLatticeFraction(identifier=name, x=a, astrand="+") 
#' barplotLatticeCounts(identifier=name,  x=a, astrand="+") 
#' 
#' @export barplotLatticeFraction
#' @export barplotLatticeCounts
#' 
NULL

#' @rdname barplot-lattice-support
barplotLatticeFraction <- function(identifier, ...) {

    if (length(list(...)) == 0) {
        e <- new.env(hash = TRUE)
    } else {
        e <- list2env(list(...))
    }

	e$ids <- unlist(e$ids)
	e$strand <- e$astrand

    if (!exists("mainvec", envir = e, inherits = FALSE)) {
		e$mainvec <- rep("",nrow(e$x))
	}
	if(class(e$mainvec)=="list"){
		e$mainvec <- unlist(e$mainvec)
	}
    if (!exists("cex.mainvec", envir = e, inherits = FALSE)) {
		e$cex.mainvec <- 1
	}
    if (!exists("main", envir = e, inherits = FALSE)) {
		e$main <- e$mainvec[e$ids %in% identifier]
	}

    if (!exists("deAnnoPlot", envir = e, inherits = FALSE)) {
        e$deAnnoPlot <- FALSE
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
    if (!exists("top.fraction.criteria", envir = e, inherits = FALSE)) {
        e$top.fraction.criteria <- "maxcount"
    }

	df <- fractionPlotDf(e$x, identifier, strand=e$strand, top.fraction.criteria=e$top.fraction.criteria)
    
    # Grah params
    my_cols <- c("green", "red")
    
    # set default values
    parset <- list(superpose.polygon=list(col=my_cols))

    scales = list(rot = c(90, 0))
    
	if (e$deAnnoPlot) {
		parset <- list(
				   layout.widths = list(
						left.padding = 0,
						axis.left = 0,
						ylab.axis.padding = 0, 
						right.padding = 0, 
						axis.right = 0
						),
				   layout.heights = list(
						top.padding = 0.1,
						between = 0.1,
						xlab.top= 0.1,
						axis.top = 0,
						main=1.1,
						main.key.padding=1,
						axis.xlab.padding = 1, 
						bottom.padding = 1, 
						axis.bottom = 0.3
						),
				   superpose.polygon=list(col=my_cols)
				   )
	
		scales = list(y = list(at = NULL, labels = NULL),
					  x = list(labels = rep("",ncol(e$x)))
					  )
					  #, rot = c(90, 0))
	}
	if(e$top.fraction.criteria=="phase"){
		df$groups <- df$phase
	}else{
		df$groups <- df$alleles
	}
    
	if(!e$middleLine) {
		b <- barchart(values ~ sample, group = groups, data = df, origin = 0, 
		    stack = TRUE, scales = scales, 
			main = list(label=e$main, cex=e$cex.mainvec), 
			ylab = e$ylab, xlab = e$xlab, 
		    par.settings = parset)
	}else if (e$middleLine) {
		b <- barchart(values ~ sample, group = groups, data = df, origin = 0, 
			stack = TRUE, scales = scales, 
			main = list(label=e$main, cex=e$cex.mainvec), 
			ylab = e$ylab, xlab = e$xlab, 
			par.settings = parset, panel=function(x, y, ...) {
				 panel.barchart(x, y, ...)
				 panel.abline(h=0.5, lty=1)
			}  )
	}else {
		stop("middleLine has to be TRUE or FALSE")
	} 
    b
    
}

#' @rdname barplot-lattice-support
barplotLatticeCounts <- function(identifier, ...) {
   
    if (length(list(...)) == 0) {
        e <- new.env(hash = TRUE)
    } else {
        e <- list2env(list(...))
    }
	
	#use requested strand instead of Gviz ranges strand
	e$strand <- e$astrand
	e$ids <- unlist(e$ids)
    #e$auto.key <- list(points = FALSE, rectangles = TRUE, space = "top", size = 2, 
	#			cex = 0.8)


    if (!exists("mainvec", envir = e, inherits = FALSE)) {
		e$mainvec <- rep("",nrow(e$x))
	}
	if(class(e$mainvec)=="list"){
		e$mainvec <- unlist(e$mainvec)
	}
    if (!exists("main", envir = e, inherits = FALSE)) {
		e$main <- e$mainvec[e$ids %in% identifier]
		#print(identifier)
		#print(e$ids)
		#print(e$main)
		#print(e$mainvec[which(e$ids %in% identifier)])
		#print(e$mainvec)
		#print(which(e$ids %in% identifier))
	}

    if (!exists("deAnnoPlot", envir = e, inherits = FALSE)) {
        e$deAnnoPlot <- FALSE
    }
    
    if (!exists("ylab", envir = e, inherits = FALSE)) {
        e$ylab <- ""
    }
    if (!exists("xlab", envir = e, inherits = FALSE)) {
        e$xlab <- ""
    }

	### Grah params set default values
	parset <- list()
	scales = list(rot = c(90, 0))
		
	makePlotDf <- function(strand){
		acounts<-  alleleCounts(e$x, strand = strand)
		arank<-  arank(e$x, strand = "*")
		#afraction<-  fraction(e$x, strand = strand)
		#amainVec<-  e$amainVec

		# prepare data to be plotted
		#a.m <- amainVec[identifier]
		a.r <- arank[[identifier]][1:2]
		a.c <- acounts[[identifier]][, a.r, drop = FALSE]
		
		values <- as.vector(t(a.c))
		sample <- vector()
		for (i in 1:nrow(a.c)) {
			sample <- c(sample, rownames(a.c)[i], rownames(a.c)[i])
		}

		if(strand=="-" && e$strand=="both"){
			values <- -values
		}
		df <- data.frame(values = values, sample = sample, alleles = rep(colnames(a.c), nrow(a.c))
)
		
		# to get right order in barchart
		df$sample <- factor(df$sample, levels = unique(df$sample))
	
		df
	}

	if(e$strand %in% c("+","-","*")){

		df <- makePlotDf(strand=e$strand)

		if (e$deAnnoPlot) {
			e$auto.key <- FALSE
			parset <- list(
					   layout.widths = list(
							left.padding = 0,
							axis.left = 0,
							ylab.axis.padding = 0, 
							right.padding = 0, 
							axis.right = 0
							),
					   layout.heights = list(
							top.padding = 0.1,
							between = 0.1,
							xlab.top= 0.1,
							axis.top = 0,
							main=1.1,
							main.key.padding=1,
							axis.xlab.padding = 1, 
							bottom.padding = 1, 
							axis.bottom = 0.3
							)
					   )
			scales = list(y = list(at = NULL, labels = NULL), rot = c(90, 0), auto.key = e$auto.key)
		}

		b <- barchart(values ~ sample, horiz = FALSE, origin = 0, group = alleles, data = df, 
			stack = FALSE, scales = scales, ylab = e$ylab, xlab = e$xlab, 
			box.ratio = 2, abbreviate = TRUE, par.settings = parset, main = e$main )

	}else if(e$strand=="both"){
		
		df <- rbind(makePlotDf(strand="+"),makePlotDf(strand="-"))	

		if (e$deAnnoPlot) {
			e$auto.key <- FALSE
			parset <- list(
					   layout.widths = list(
							left.padding = 0,
							axis.left = 0,
							ylab.axis.padding = 0, 
							right.padding = 0, 
							axis.right = 0
							),
					   layout.heights = list(
							top.padding = 0.1,
							between = 0.1,
							xlab.top= 0.1,
							axis.top = 0,
							main=1.1,
							main.key.padding=1,
							axis.xlab.padding = 1, 
							bottom.padding = 1, 
							axis.bottom = 0.3
							)
					   )
			scales = list(y = list(at = NULL, labels = NULL), rot = c(90, 0), auto.key = e$auto.key)
		}
		
		b <- barchart(values ~ sample, horiz = FALSE, origin = 0, group = alleles, data = df, 
			stack = FALSE, scales = scales, ylab = e$ylab, xlab = e$xlab, 
			box.ratio = 2, abbreviate = TRUE, par.settings = parset, main = e$main)
	}else {
		stop("strand must be + - * or both")
	}
    b
}



#' coverage matrix of GAlignmentsList
#' 
#' Get coverage per nucleotide for reads covering a region
#' 
#' a convenience function to get the coverage from a list of reads stored in
#' GAlignmnetsList, and returns by default a list with one matrix, and
#' information about the genomic start and stop positions.
#' 
#' @param BamList GAlignmentsList containing reads over the region to calculate
#' coverage
#' @param strand strand has to be '+' or '-'
#' @param ignore.empty.bam.row argument not in use atm
#' @author Jesper R. Gadin
#' @keywords coverage
#' @examples
#' 
#' r <- reads
#' seqlevels(r) <- '17'
#' covMatList <- coverageMatrixListFromGAL(BamList=r, strand='+')
#' 
#' @export coverageMatrixListFromGAL
coverageMatrixListFromGAL <- function(BamList, strand = "*", ignore.empty.bam.row = TRUE) {
    
    # If having common start and end points for all gviz track objects the matrix
    # will start on the specific start regardless if there are reads in the bamList
    # or not.
    
    # TODO, to conveniently access data without loading into memory, the bam file
    # should be read again by using the argument BamPath.
    
    GAL <- BamList
    
    # Could be good with a check that the matrix is not longer than CNTNAP2 -
    # 2300000bp long which is the longest gene But will wait with that.
    pstrand = FALSE
    mstrand = FALSE

    if (!is.null(strand)) {
        if (strand == "+") {
            pstrand = TRUE
        } else if (strand == "-") {
            mstrand = TRUE
        } else if (strand == "*") {
            mstrand = TRUE
            pstrand = TRUE
        } else if (strand == "both") {
            mstrand = TRUE
            pstrand = TRUE
        } else {
            stop("strand has to be '+' or '-' if not NULL\n")
        }
    }
    
    if (!length(seqlevels(GAL)) == 1) {
        stop("can only be one seq level\n")
    }
    
    # get start and end before filtering on strand, will make things easier
    # downstream.
    suppressWarnings(bamStart <- min(min(start(GAL))))
    suppressWarnings(bamEnd <- max(max(end(GAL))))
    bamWidth <- bamEnd - bamStart + 1
    
    # if(is.null(start) | is.null(end)){
    start <- bamStart
    end <- bamEnd
    width <- bamWidth
    # }
    
    if (pstrand) {
        GALp <- GAL[strand(GAL) == "+"]
    }
    if (mstrand) {
        GALm <- GAL[strand(GAL) == "-"]
    }
    
    if (pstrand) {
        matP <- matrix(0, ncol = (width), nrow = length(GAL))
    }
    if (mstrand) {
        matM <- matrix(0, ncol = (width), nrow = length(GAL))
    }
    if (pstrand) {
        rownames(matP) <- names(GAL)
    }
    if (mstrand) {
        rownames(matM) <- names(GAL)
    }
    ################################################## 
    
    covVecFromGA <- function(GA) {
        mcols(GA) <- NULL
        one <- unlist(grglist(GA))
        covRle <- coverage(one)[[1]]
        cov <- as.integer(window(covRle, start, end))
        cov
    }
    
    if (pstrand) {
        for (i in 1:length(GALp)) {
            matP[i, ] <- covVecFromGA(GALp[[i]])
        }
    }
    if (mstrand) {
        for (i in 1:length(GALm)) {
            matM[i, ] <- covVecFromGA(GALm[[i]])
        }
    }
    
    # make mat from matP or matM
    if (strand=="+") {
        mat <- matP
    }
    if (strand=="-") {
        mat <- matM
    }
    if (strand=="*") {
        mat <- matM+matP
    }
    
    if (strand=="both") {
        mat <- list(plus=matP,minus=matM)
    }

    # store in a list
    if (!is.null(strand)) {
        retList <- list(mat, start, end)
    } else {
        stop("strand must be present")
    }
    # set name on list
    if (!is.null(strand)) {
        names(retList) <- c("mat", "start", "end")
    } else {
        stop("strand must be present")
    }
    
    retList
}



#' implode list of arguments into environment
#' 
#' apply on list of variables to be put in the local environment
#' 
#' help the propagation of e.g. graphical paramters 
#' 
#' @param x list of variables
#' @author Jesper R. Gadin
#' @keywords implode
#' @examples
#' 
#' lst <- list(hungry='yes', thirsty='no')
#' implodeList(lst)
#' #the check ls()
#'  ls()
#' @export implodeList
implodeList <- function(x) {
    oname <- deparse(substitute(x))
    eval(parse(text = paste0("for(i in 1:length(", oname, ")){assign(names(", oname, 
        ")[i],", oname, "[[i]])}")), parent.frame())
} 

#' alleleCounts from bam file
#' 
#' count alleles before creating ASEse.
#' 
#' counts the alleles in a bam file based on GRanges positions. 
#' 
#' 
#' @param gr GRanges that contains SNPs of interest
#' @param pathToDir path to directory of bam files
#' @param flag specify one flag to use as filter, default is no filtering. 
#' allowed flags are 99, 147, 83 and 163
#' @param scanBamFlag set a custom flag to use as filter
#' @param return.class type of class for the returned object
#' @param verbose makes funciton more talkative
#' @author Jesper R. Gadin
#' @keywords allelecount counting
#' @examples
#'
#' data(GRvariants)
#' gr <- GRvariants
#'
#' ##not run at the moment
#' #pathToDir <- system.file('inst/extdata/ERP000101_subset', package='AllelicImbalance')
#' #ar <- countAllelesFromBam(gr, pathToDir)
#'  
#' @export countAllelesFromBam
countAllelesFromBam <- function(gr, pathToDir, flag=NULL, scanBamFlag=NULL, return.class="array", verbose=TRUE) {

	bamDir <- normalizePath(pathToDir)
	allFiles <- list.files(bamDir, full.names = TRUE)
	bamFiles <- allFiles[grep(".bam$", allFiles)]
	if (length(bamFiles) == 0) {
		stop(paste("No bam files found in", bamDir))
	}
	if (!all(file.exists(paste(bamFiles, ".bai", sep = "")))) {
		if (verbose) {
			cat(paste("The bam files in UserDir are required to also have", ".bam.bai index files.", 
				" Trying to run indexBam function on each", "\n"), )
		}
		indexBam(bamFiles)
		if (!all(file.exists(paste(bamFiles, ".bai", sep = "")))) {
			stop("The bam files in UserDir are required to also have", ".bam.bai index files.")
		} else {
			if (verbose) {
				cat(paste("Succesfully indexed all bamFiles in UserDir", UserDir, 
				  "\n"))
			}
		}
	}
	
	#scanBamFlag
	if(is.null(scanBamFlag) & is.null(flag)){
		flag <- scanBamFlag()
	}
	if(!is.null(scanBamFlag) ){
		if(length(scanBamFlag)==2){
			flag <- scanBamFlag
		}else{
			stop("scanBamFlag has to be the return values from scanBamFlag()")
		}
	}
	#flags can be 99 147 83 or 163
	if(length(flag)==1 & is.null(scanBamFlag)){
		if(! (flag%in%c(99,147,83,163))){
			stop(paste("flag values can only be 99 147 83 or 163",
					   "the input flag value was", flag,sep=" "))
		}
		
		if(flag==99){

			flag=scanBamFlag( isPaired = TRUE,
					isProperPair = TRUE,
					isFirstMateRead = TRUE,	
					isMateMinusStrand = TRUE
					)
		}else if(flag==83){

			flag=scanBamFlag(
					isPaired = TRUE,
					isProperPair = TRUE,
					isFirstMateRead = TRUE,	
					isMinusStrand = TRUE
					)
		}else if(flag==147){
		
			flag=scanBamFlag(
					isPaired = TRUE,
					isProperPair = TRUE,
					isFirstMateRead = FALSE,	
					isMinusStrand = TRUE
					)
		}else if(flag==163){
		
			flag=scanBamFlag(
					isPaired = TRUE,
					isProperPair = TRUE,
					isFirstMateRead = FALSE,	
					isMateMinusStrand = TRUE
					)
		}
	}


	if(verbose){
		cat("sam flag used:",flag,"\n")
	}
	fls <- PileupFiles(bamFiles)
	
	countF <-
		function(x){
		x[["seq"]][-5,,1]

	}                     

	which <- gr
	p1 <- ApplyPileupsParam(flag=flag,
							which=which, 
							minBaseQuality = 0L,
							what="seq",
							yieldBy = "position",
							yieldAll=TRUE
							)

	res <- applyPileups(fls, countF, param=p1)

	ar <- array(unlist(res), dim=c(4, length(fls), length(which)),
				dimnames=list(c("A","C","G","T"), 
							  names(fls), 
							  names(which)))
	ar <- aperm(ar,dim=c(3,2,1))
	
	ar
}

#' ASEset from bam file
#' 
#' count alleles and create an ASEset direct from bam file instead of reading into R first.
#' 
#' counts the alleles in a bam file based on GRanges positions. 
#' 
#' 
#' @param gr GenomicRanges of SNPs to create ASEset for
#' @param PE if paired end or not (default: TRUE)
#' @param pathToDir Directory of bam files with index in same directory
#' @param strandUnknown default: FALSE
#' @param ... passed on to countAllelesFromBam function
#' @param flagsMinusStrand flags that mark reads coming from minus strand
#' @param flagsPlusStrand flags that mark reads coming from plus strand
#' @author Jesper R. Gadin
#' @keywords ASEset
#' @examples
#'
#' data(GRvariants)
#' gr <- GRvariants
#'
#' ##no execution at the moment
#' #pathToDir <- system.file('inst/extdata/ERP000101_subset', package='AllelicImbalance')
#' #a <- ASEsetFromBam(gr, pathToDir)
#'  
#' @export ASEsetFromBam

ASEsetFromBam <- function(gr, pathToDir,PE=TRUE, flagsMinusStrand=c(83,163), flagsPlusStrand=c(99,147), strandUnknown=FALSE, ...) {

	if(!PE){
		stop("no support for SE atm")
	}

	if(PE==TRUE){
		#minus strand
		arm1 <- countAllelesFromBam(gr, pathToDir, flag=83)
		arm2 <- countAllelesFromBam(gr, pathToDir, flag=163)
		arm <- arm1 + arm2

		#plus strand
		arp1 <- countAllelesFromBam(gr, pathToDir, flag=99)
		arp2 <- countAllelesFromBam(gr, pathToDir, flag=147)
		arp <- arp1 + arp2
	}
	
	#ASEsetFromArray
	if(!strandUnknown){
		a <- ASEsetFromArrays(gr, countsPlus = arp, 
			countsMinus = arm)

	}else{	
		a <- ASEsetFromArrays(gr, countsUnknown = arp+arm) 
	}
	a
}


#' makes masked fasta reference
#' 
#' Replaces all selected positions in a fasta file with the character N
#' 
#' @param fastaIn character string of the path for the fasta file to be used
#' @param fastaOut character string of the path for the masked fasta file (no extension)
#' @param posToReplace GRanges object with the genomic ranges to replace 
#' @param splitOnSeqlevels write on file for each seqlevel to save memory
#' @param verbose makes function more talkative
#' @author Jesper R. Gadin
#' @keywords masked fasta reference
#' @examples
#'
#' data(ASEset.sim)
#' gr <- rowRanges(ASEset.sim) 
#' fastaIn <- system.file('extdata/hg19.chr17.subset.fa', package='AllelicImbalance')
#' makeMaskedFasta(fastaIn=fastaIn, fastaOut="fastaOut",posToReplace=gr) 
#' 
#'  
#' @export makeMaskedFasta
makeMaskedFasta <- function(fastaIn, fastaOut, posToReplace, splitOnSeqlevels=TRUE, verbose=TRUE){

	#does the inFasta file exist?
	if(!file.exists(fastaIn)){
		stop("fasta infile doesnt exist")
	}

	#make FaFile
	fl <- FaFile(fastaIn)

	#check if the index file is present otherwise tell the user to use the
	#indexFa(FaFile("pathToReference")) command
	if(!file.exists(paste(fastaIn,".fai",sep=""))){
		cat("could not find index file\n")
		cat("creates a new index file")
		indexFa(FaFile(fastaIn))
		cat("finished creating new index file")
	}
	#indexFa(fl) #creates a new file as index in the same directory but
	#with extension *.fai 

	#open,scan,close file
	open(fl)
	#check if seqleveles are present
	fa.info <- scanFaIndex(fl)
	if(!(all(seqlevels(posToReplace) %in% seqlevels(fa.info)))){
		close(fl)
		stop("seqlevels in object x are not in fasta index file")
	}

	###############LOOOP
	chrs <- seqlevels(posToReplace)
	for (chr in chrs){
	#	chr <- "chr1"

		searchArea <- fa.info[seqnames(fa.info)==chr]
		seq <- scanFa(fl,param=searchArea)

		#replace the SNPs with N
		toReplace <- unique(as.integer(ranges(posToReplace)))
		seq <- seq[[1]]
		seq[toReplace] <- "N"

		if(verbose){cat("replaced", length(toReplace),"instances with N\n")}

		#write new file
		#library("seqinr")
		outfile <- paste(fastaOut,chr,".fa",sep="")

		write.fasta(as.character(seq),names=chr, file.out=outfile, nbchar=80)	
		if(verbose){cat("wrote chr",chr,"to file\n")}

		gc()
	}
	close(fl)
	if(verbose){cat("all chromosomes written to file\n")}
}



#' global analysis wrapper
#' 
#' A wrapper to make a global analysis based on paths for BAM, VCF and GFF files
#' 
#' @param pathBam path to bam file
#' @param pathVcf path to vcf file
#' @param pathGFF path to gff file
#' @param verbose makes function more talkative
#' @author Jesper R. Gadin
#' @keywords global wrapper
#' @examples
#'
#' #empty as function doesn't exist
#' 
#' @export 
gba <- function(pathBam,pathVcf,pathGFF=NULL, verbose){

	#summarize counts
	
	#detectAI

	
}



