#get SnpAlleleCounts
setMethod("alleleCounts",
    signature(x = "ASEset"),
    function (x,strand="nonStranded")
	{
		if(!sum(strand %in% c("+","-","*","nonStranded"))>0){stop("strand parameter has to be either '+', '-', '*' or 'nonStranded' ")}
	
		if(strand=="+"){
			el <- "countsPlus"
		}else if(strand=="-"){
			el <- "countsMinus"
		}else if(strand=="*"){
			el <- "countsUnknown"
		}else if(strand=="nonStranded"){
			el <- "countsNonStranded"
		}else{stop("unknown strand option")}	

		#check if strand option is present as assay
		if(!(el%in%names(assays(x)))){stop("strand is not present as assay in ASEset object")}

		#assume alleleCount information is
		#stored as element 1
		alleleCountList <- list()

		for(i in 1:nrow(assays(x)[[el]])){
			mat <- assays(x)[[el]][i,,]
			if(class(mat)=="numeric"){
					mat <- t(mat)
					colnames(mat) <- x@variants
				}else{
					colnames(mat) <- x@variants
				}
			rownames(mat) <- colnames(x)
			alleleCountList[[i]] <- mat
		}	
		#add snp id
		names(alleleCountList) <- rownames(x)	
	
		#return object
		alleleCountList

	}
)

setMethod("mapBias",
    signature(x = "ASEset"),
    function (x)
	{
		#assume alleleCount information is
		#stored as element 1
		mapBiasList <- list()
		for(i in 1:nrow(x)){
			mat <- assays(x)[["mapBias"]][i,,]
			if(class(mat)=="numeric"){
				dim(mat) <- c(1,5)
				rownames(mat) <- colnames(x)
			}
			colnames(mat) <- x@variants

			mapBiasList[[i]] <- mat
		}	
		#add snp id
		names(mapBiasList) <- rownames(x)	
		

		#return object
		mapBiasList
    
	}
)

setMethod("fraction",
    signature(x = "ASEset"),
    function (x, strand="nonStranded",verbose=FALSE)
	{

		if(!sum(strand %in% c("+","-","*","nonStranded"))>0){stop("strand parameter has to be either '+', '-', '*' or 'nonStranded' ")}
	
		if(strand=="+"){
			el <- "countsPlus"
		}else if(strand=="-"){
			el <- "countsMinus"
		}else if(strand=="*"){
			el <- "countsUnknown"
		}else if(strand=="nonStranded"){
			el <- "countsNonStranded"
		}else{stop("unknown strand option")}	

		#check if strand option is present as assay
		if(!(el%in%names(assays(x)))){stop("strand is not present as assay in ASEset object")}
		
		fractionList <- list()

		for(i in 1:nrow(x)){
			#getting revelant data
			tmp <- alleleCounts(x, strand)[[i]]
			
			#calculating major and minor allele, warn if the two remaining alleles have too many counts
			if(nrow(tmp)>1){
				countsByAllele<-apply(tmp,2,sum,na.rm=TRUE)
			}else{
				countsByAllele<-tmp
			}
			majorAllele <-colnames(tmp)[order(countsByAllele,decreasing=TRUE)][1]
			minorAllele <-colnames(tmp)[order(countsByAllele,decreasing=TRUE)][2]
			majorAndMinorFraction <- sum(countsByAllele[c(majorAllele,minorAllele)]) / sum(countsByAllele)
			if(verbose & majorAndMinorFraction < 0.9){
				cat(paste("Snp","was possible tri-allelic, but only two most frequent alleles were plotted. Counts:"),"\n")
				cat(paste(paste(names(countsByAllele),countsByAllele,sep="="),collapse=", "),"\n")
			}
			
			#calculating percentage and ylim and setting no-count samples to colour grey90
			fraction <- tmp[,majorAllele] / (tmp[,majorAllele] + tmp[,minorAllele])
			fraction[is.nan(fraction)]<-1
			fractionList[[i]] <- fraction
		}
		names(fractionList) <- rownames(x)	
		
		#return object
		#fractionList
		as.matrix(as.data.frame(fractionList))
    
	}
)




