#
#.boxplot.internal <- function(x, strand="*",
#	labels.axis=TRUE, ...)
#{
#
#	bp <- boxplot(frequency(x,strand=strand), ...)
#	
#	#add numbers of used SNPs ontop of boxplot
#	labels <- hetPerSample(x,strand=strand)
#	axis(side=3,labels=labels,at=1:length(labels))
#
#	#add axis labels
#	if(labels.axis){
#		mtext("Nr. Of Heterozygots Used",side=3, padj=-3, cex=1, font=2)
#		mtext("Reference Allele Frequency", side=2,padj=-3, cex=1, font=2)
#		mtext("Samples",side=1,padj=3, cex=1, font=2)
#	}
#
#	#add red line for 0.5
#	abline(h=0.5, col="red")
#	
#	invisible(bp)
#}
#
#
#
