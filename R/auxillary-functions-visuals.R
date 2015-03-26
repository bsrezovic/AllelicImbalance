#' add legend to AllelicImbalance barplot
#' 
#' adds a very customizable legend function for AllelicImbalance barplots.
#' 
#' the function is preferably called from within the AllelicImbalance barplot method.
#' 
#' 
#' @param lowerLeftCorner position of the plot to add legend to (default c(0,0))
#' @param size scale the plot, default is 1
#' @param rownames rownames in legend
#' @param colnames colnames in legend
#' @param boxsize size of each box fill
#' @param boxspace space inbetween the box fill
#' @param fgCol color for allele1
#' @param bgCol color for allele2
#' @param ylegendPos placement of the legend within the plot for y
#' @param xlegendPos placement of the legend within the plot for x
#' @param cex size of legend text
#' @author Jesper R. Gadin
#' @keywords legend barplot
#' @examples
#'
#' #code placeholders
#' #< create a barplot with legend >
#' #< add legend >
#'  
#' @export legendBarplot

legendBarplot <- function(lowerLeftCorner, size, rownames, colnames, boxsize=1, boxspace=1,fgCol,bgCol,
						  ylegendPos=1, xlegendPos=0.96, cex=1){
	cex.legend <- cex

	#first fill/box
	x = lowerLeftCorner[1] + size[1] * seq(xlegendPos, (xlegendPos - 
		(0.02 * boxspace * (length(unique(fgCol)) - 1))),
			length = length(unique(fgCol)))
	y = (lowerLeftCorner[2] + size[2] * (ylegendPos + 0.02 * boxspace)) * rep(1,
			length(unique(fgCol)))
	squares	= rep(c(size[1] * 0.012 * boxsize), length(unique(fgCol)))

		symbols(x=x, y = y , bg = unique(fgCol), squares = squares, add = TRUE, 
			inches = FALSE,xpd=TRUE)

	if(length(rownames)==2){
		#second fill/box
		x = lowerLeftCorner[1] + size[1] * seq(xlegendPos, (xlegendPos - 
			(0.02 * boxspace * (length(unique(bgCol)) - 1))),
				length = length(unique(bgCol)))
		y=(lowerLeftCorner[2] + size[2] * 1) * rep(ylegendPos, length(unique(bgCol)))
		squares=rep(c(size[1] * 0.012 * boxsize), length(unique(bgCol)))

		symbols(x=x, y = y , bg = unique(bgCol), squares = squares, add = TRUE, 
			inches = FALSE,xpd=TRUE)
	}
	# row-lab
	if(length(rownames)==2){

		x = c(lowerLeftCorner[1] + (size[1] * (c(xlegendPos, xlegendPos)+0.02 * boxspace)))
		y = lowerLeftCorner[2] + (size[2] * c(ylegendPos+(0.02 * boxspace), 1))
		text(x=x, 
		  y = y, rownames, 
		  srt = 0, cex = cex.legend, adj = c(0, 0.5), xpd = TRUE)

	}else if(length(rownames)==1){

		x = c(lowerLeftCorner[1] + (size[1] * (c(xlegendPos)+0.02 * boxspace)))
		y = lowerLeftCorner[2] + (size[2] * c(ylegendPos+(0.02 * boxspace)))
		text(x=x, 
		  y = y, rownames, 
		  srt = 0, cex = cex.legend, adj = c(0, 0.5), xpd = TRUE)

	}
	# col-lab
	  x = lowerLeftCorner[1] + size[1] * seq(xlegendPos, (xlegendPos - ((0.02 * boxspace) * 
	  (length(unique(fgCol)) - 1))), length = length(unique(fgCol)))
	  y = lowerLeftCorner[2] + size[2] * (ylegendPos+0.03*boxspace)
	text(x=x, 
		y = y, colnames, 
		srt = 90, cex = cex.legend, adj = c(0, 0.5), xpd = TRUE)
}

#' add annotation to AllelicImbalance barplot
#' 
#' adds a customizable annotation functionality for AllelicImbalance barplots.
#' 
#' the function is preferably called from within the AllelicImbalance barplot method.
#' 
#' 
#' @param strand strand, "+", "-", "*" or "both"
#' @param snp integer for the described snp
#' @param lowerLeftCorner position of the plot to add legend to (default c(0,0))
#' @param annDfPlus annotation dataframe plus strand
#' @param annDfMinus annotation dataframe minus strand
#' @param cex size of annotation text
#' @param ypos relative y-axis position for the annotation text 
#' @param interspace space between each annotation block
#' @author Jesper R. Gadin
#' @keywords annotation barplot
#' @examples
#'
#' #code placeholders
#' #< create a barplot without annotation >
#' #< add annotation >
#'  
#' @export annotationBarplot

annotationBarplot <- function(strand,snp, lowerLeftCorner, annDfPlus, annDfMinus,
							  cex=0.7, ypos=0, interspace=1){
	i <- snp

	#ypos, is relative to the default

	  # check which columns to extract
	  dfType <- c("symbol", "exon_id", "tx_id", "cds_id")
	  TFann <- (dfType %in% colnames(annDfPlus)) | (dfType %in% colnames(annDfMinus))


	if (strand == "+" | strand == "-" | strand == "*" ) {

	  # placed on left side
	  if (strand == "+" | strand == "*" ) {
		  df <- annDfPlus[, TFann, drop = FALSE]
		
		yPosTmp <- lowerLeftCorner[2] + (1.03) + ypos
		text <- paste("Plus-strand")
		text(x = lowerLeftCorner[1] + 0.01, y = yPosTmp, labels = text, 
		  adj = 0, xpd = TRUE, cex=cex)
		yPosTmp <- yPosTmp - (0.03*interspace)
		
		for (j in 1:length(colnames(df))) {
		  text <- paste(colnames(df)[j], df[i, j], sep = ": ")
		  text(x = lowerLeftCorner[1] + 0.01, y = yPosTmp, labels = text, 
			adj = 0, xpd = TRUE, cex=cex)
		  yPosTmp <- yPosTmp - (0.02*interspace)
		}
	  }
	  # placed on left side
	  if (strand == "-" & (!strand == "*" )) {
	  df <- annDfMinus[, TFann, drop = FALSE]
		
		yPosTmp <- lowerLeftCorner[2] + (1.03) + ypos
		text <- paste("Minus-strand")
		text(x = lowerLeftCorner[1] + 0.01, y = yPosTmp, labels = text, 
		  adj = 0, xpd = TRUE, cex=cex)
		yPosTmp <- yPosTmp - (0.03*interspace)
		for (j in 1:length(colnames(df))) {
		  text <- paste(colnames(df)[j], df[i, j], sep = ": ")
		  text(x = lowerLeftCorner[1] + 0.01, y = yPosTmp, labels = text, 
			 adj = 0, xpd = TRUE, cex=cex)
		  yPosTmp <- yPosTmp - (0.02*interspace)
		}
	  }
	  # placed on right side
	  if (strand == "*" ) {
	  df <- annDfMinus[, TFann, drop = FALSE]
		
		yPosTmp <- lowerLeftCorner[2] + (1.03) + ypos
		text <- paste("Minus-strand")
		text(x = lowerLeftCorner[1] + 1, y = yPosTmp, labels = text, 
		  adj = 0.95, xpd = TRUE, cex=cex)
		yPosTmp <- yPosTmp - (0.03*interspace)
		
		for (j in 1:length(colnames(df))) {
		  text <- paste(df[i, j], colnames(df)[j], sep = " :")
		  text(x = lowerLeftCorner[1] + 1, y = yPosTmp, labels = text, 
			 adj = 0.95, xpd = TRUE, cex=cex)
		  yPosTmp <- yPosTmp - (0.02*interspace)
		}
	  }
	

	########
	# Dual barplot 
	########
	}else if (strand == "both") {

	  
	  # plus strand (over)
	  df <- annDfPlus[, TFann, drop = FALSE]
	  yPosTmp <- lowerLeftCorner[2] + (1.01) + ypos
	  for (j in 1:length(colnames(df))) {
		text <- paste(colnames(df)[j], df[i, j], sep = ": ")
		text(x = lowerLeftCorner[1] + 0.01, y = yPosTmp, labels = text, 
		   adj = 0, cex=cex)
		yPosTmp <- yPosTmp - (0.02*interspace)
	  }
	  
	  # minus strand (under)
	  df <- annDfMinus[, TFann, drop = FALSE]
	  
	  yPosTmp <- lowerLeftCorner[2] - (0.01) + ypos
	  for (j in 1:length(colnames(df))) {
		text <- paste(colnames(df)[j], df[i, j], sep = ": ")
		text(x = lowerLeftCorner[1] + 0.01, y = yPosTmp, labels = text, 
		  adj = 0, cex=cex)
		yPosTmp <- yPosTmp + (0.02*interspace)
	  }
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


