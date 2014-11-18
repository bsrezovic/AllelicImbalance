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

