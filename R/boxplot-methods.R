#' boxplots 
#' 
#' uses base graphics box plot
#' 
#' The boxplot will show the density over frequencies for each sample
#' 
#' @name boxplot
#' @rdname boxplot
#' @aliases boxplot boxplot,ReferenceBias-method 
#' @docType methods
#' @param x \code{ReferenceBias} object
#' @param strand '+','-' or '*'
#' @param labels.axis logical
#' @param ... arguments to forward to internal boxplots function
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords plot box
#' @examples
#' 
#' #load example data
#' 
#' data("ReferenceBias")
#' boxplot(ReferenceBias)
#' 
NULL

#' @rdname boxplot
setMethod("boxplot", signature(x = "ReferenceBias"), function(x, strand="*",
	labels.axis=TRUE, ...){
	bp <- boxplot(frequency(x,strand=strand), ...)
	
	#add numbers of used SNPs ontop of boxplot
	labels <- hetPerSample(x,strand=strand)
	axis(side=3,labels=labels,at=1:length(labels))

	#add axis labels
	if(labels.axis){
		mtext("Nr. Of Heterozygots Used",side=3, padj=-3, cex=1.5, font=2)
		mtext("Reference Allele Frequency", side=2,padj=-3, cex=1.5, font=2)
		mtext("Samples",side=1,padj=3, cex=1.5, font=2)
	}

	#add red line for 0.5
	abline(h=0.5, col="red")
	
	invisible(bp)
})



