#setMethods("initialize") is not required at the moment

#wrapper function to create ASEset object from AlleleCountList
ASEset <- function(sset,variants){
	
	#create object
	new("ASEset",sset,variants=variants)
}	

ASEsetFromCountList <- function(rowData, countListNonStranded=NULL, countListPlus=NULL, countListMinus=NULL, countListUnknown=NULL, colData=NULL, mapBiasExpMean=NULL,verbose=FALSE ,...){	
	
	if(verbose){
		cat("rowData\n")
		cat(class(rowData))
		cat("countListNonStranded\n")
		cat(class(countListNonStranded))
		cat("countListPlus\n")
		cat(class(countListPlus))
		cat("countListMinus\n")
		cat(class(countListMinus))
		cat("countListUnknown\n")
		cat(class(countListUnknown))
		cat("colData\n")
		cat(class(colData))
		cat("mapBiasExpMean\n")
		cat(class(mapBiasExpMean))
	}
	#check that at least one of the countList options are not null
	if(is.null(c(countListNonStranded,countListPlus,countListMinus,countListUnknown))){
		stop("at least one of the countList options has to be specified")
	}

	countLists <- c("countListNonStranded","countListPlus","countListMinus","countListUnknown")[(c(!is.null(countListNonStranded),!is.null(countListPlus),!is.null(countListMinus),!is.null(countListUnknown)))] 		
	
	#check that all lengths are the same in all lists
	l <- unlist(lapply(countLists,function(x){length(get(x))}))
	if(sum(!(l == l[1]))>0){stop("one or more list contains more or less list elements than the others")}
	
	#check that the number of columns in each dataframe are the same
	l <- unlist(lapply(countLists,function(x){lapply(get(x),ncol)}))
	m <- matrix(unlist(l),ncol=length(l))
	if(sum(!(m == m[1]))>0){stop("one or more list contains more or less columns than the others")}
	
	#check that the number of rows in each dataframe are the same
	l <- unlist(lapply(countLists,function(x){lapply(get(x),nrow)}))
	m <- matrix(unlist(l),ncol=length(l))
	if(sum(!(m == m[1]))>0){stop("one or more list contains more or less rows than the others")}

	#check that all colnames in all dataframes are the same in same order
	cMatrix <- matrix(NA,nrow=length(unlist(unique(lapply(countLists[1], colnames)))),ncol=length(countLists))
	for( i in 1:length(countLists)) { #within list comparision
		countListName <- countLists[i]
		countList <- get(countListName)
		l <- lapply(countList, colnames)
		m <- matrix(unlist(l),ncol=length(l))
		if(sum(!(m==m[,1]))>0){stop(paste("the names or order of names are not the same in all data frames in list",countListName))}
		cMatrix[,i] <- m[,1]
	}
	if(sum(!(cMatrix==cMatrix[,1]))>0){stop(paste("the names or order of names are not the same between all data frame lists"))}

	#check that all rownames in all dataframes are the same in same order
	cMatrix <- matrix(NA,nrow=length(unlist(unique(lapply(countLists[1], rownames)))),ncol=length(countLists))
	for( i in 1:length(countLists)) { #within list comparision
		countListName <- countLists[i]
		countList <- get(countListName)
		l <- lapply(countList, rownames)
		m <- matrix(unlist(l),ncol=length(l))
		if(sum(!(m==m[,1]))>0){stop(paste("the names or order of names are not the same in all data frames in list",countListName))}
		cMatrix[,i] <- m[,1]
	}
	if(sum(!(cMatrix==cMatrix[,1]))>0){stop(paste("the names or order of names are not the same between all data frame lists"))}

	#check mapBiasExpMean
	if(!(is.null(mapBiasExpMean))) {
		if(!class(mapBiasExpMean)=="array"){stop("mapBiasExpMean has to be of class array")}
		if(!length(dim(mapBiasExpMean))==3){stop("mapBiasExpMean has to have three dimensions")}
	
		for(i in 1:(dim(m)[2])){
			#check that no snp and sample exceeds sum 1 frequency
			if(!sum(!apply(m[,i,],1,sum)==1)==0){stop(paste("for each snp and sample the sum of allele frequencies must sum to one,
				for sample element",i,"and snp nr",paste(which(!apply(m[,i,],1,sum)==1),collapse=" "),"this was not fullfilled"))}
			#check for tri-allelic cases, which we dont allow as mapping biases.
			if(!sum(!(apply((m[,i,]>0),1,sum)==2))==0){stop(paste("tri-allelic SNPs have not been implemented yet. Please write an email if this is of interest.\n 
				Tri-allelic case was found for sample nr",i,"and snp nr",paste(which(!apply((m[,i,]>0),1,sum)==2),collapse=" ")))}
		}
	}



	#choose a common countList by picking the first one, for dimension info
	countList <- get(countLists[1])
	ind <- length(unlist(unique(lapply(countList, rownames))))
	snps <- length(countList)

	#SimpleList init
	assays <- SimpleList()

	#plus
	if(!is.null(countListPlus)){
		ar1 <- array(NA,c(snps,ind,5)) #empty array that handles only four nucleotides + one del columns
		for(i in 1:snps){
			ar1[i,,] <- countListPlus[[i]]
		}
		assays[["countsPlus"]] <- ar1
	}
	#minus
	if(!is.null(countListMinus)){
		ar2 <- array(NA,c(snps,ind,5)) #empty array that handles only four nucleotides + one del columns
		for(i in 1:snps){
			ar2[i,,] <- countListMinus[[i]]
		}
		assays[["countsMinus"]] <- ar2

	}
	#unknown
	if(!is.null(countListUnknown)){
		ar3 <- array(NA,c(snps,ind,5)) #empty array that handles only four nucleotides + one del columns
		for(i in 1:snps){
			ar3[i,,] <- countListUnknown[[i]]
		}
		assays[["countsUnknown"]] <- ar3

	}
	#nonStranded
	if(!is.null(countListNonStranded)){
		ar4 <- array(NA,c(snps,ind,5)) #empty array that handles only four nucleotides + one del columns
		for(i in 1:snps){
			ar4[i,,] <- countListNonStranded[[i]]
		}
		assays[["countsNonStranded"]] <- ar4
	}
	
	#assign mapBiasExpMean
	if(is.null(mapBiasExpMean)) {
		assays[["mapBias"]] <-  getDefaultMapBiasExpMean3D(countList)
	}else{
		assays[["mapBias"]] <- mapBiasExpMean
	}

	if(is.null(colData)){
		colData <- DataFrame(row.names=unlist(unique(lapply(countList, rownames))))
	}

	sset <- SummarizedExperiment(assays=assays,rowData=rowData,colData=colData,...)

	rownames(sset) <- names(countList)


	#use colnames in list matrices as variants
	variants <-unlist(unique(lapply(countList,colnames)))
	
	#create object
	ASEset(sset,variants=variants)
}


