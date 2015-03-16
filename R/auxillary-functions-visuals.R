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

