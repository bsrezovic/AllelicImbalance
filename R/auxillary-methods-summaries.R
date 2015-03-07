
#function to summarize selected variable combinations.
#x against all other variables

#generate example
#data(ASEset)
#object <- detectAI(ASEset, 
#			threshold.count.sample=1:50,
#			threshold.frequency=seq(0,0.5,by=0.01),
#			threshold.delta.frequency=seq(0,0.5,by=0.01),
#			threshold.pvalue=rev(seq(0.001,0.05, by=0.005))
#)

frequency_vs_threshold_variable_summary <- function(object,var="threshold.count.sample"){

	#check if assay is present
	fr <- assays(object)[["reference.frequency"]]
	ar.var <- assays(object)[[var]]
	ar.fr <- array(fr,dim=c(nrow(fr),ncol(fr),dim(ar.var)[3]),
				  dimnames=list(rownames(object),colnames(object),NULL) )

	is.na(ar.var) <- FALSE 
	ar.fr[!ar.var] <- NA

	#colSums(ar.fr,na.rm=TRUE)
	apply(ar.fr,c(2, 3),mean,na.rm=TRUE)
}

detectedAI_vs_threshold_variable_summary <- function(object,var="threshold.count.sample"){

	#check if assay is present

	apply(assays(object)[[var]],c(2, 3),sum,na.rm=TRUE)

}


