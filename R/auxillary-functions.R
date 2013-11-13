#These are non class specific functions, availanble as utilities to ease for example import of sequence data.

impBamGRL <- function(UserDir,searchArea,verbose=TRUE){
	#Set parameters
	which <- searchArea #A GRanges, RangesList, RangedData, or missing object, from which a IRangesList instance will be constructed.
	what <- scanBamWhat() #A character vector naming the fields to return. scanBamWhat() returns a vector of available fields. Fields are described on the scanBam help page.
	flag <- scanBamFlag(isUnmappedQuery=FALSE)	
	param<-ScanBamParam(flag=flag,which=which,what=what) #store ScanBamParam in param.
	
	#Point to correct directory and create a BamFileList object
	bamDir <- normalizePath(UserDir) #Point to the directory containing your Bam files and its respective bam.bai files.
	allFiles <- list.files(bamDir,full.names = TRUE) #list files in a folder.
	bamFiles <- allFiles[grep(".bam$",allFiles)] 	#list only the files ending in .bam .
	if(length(bamFiles) == 0)stop(paste("No bam files found in",bamDir))
	if(!all(file.exists(paste(bamFiles,".bai",sep="")))){
		if(verbose)cat(paste("The bam files in UserDir are required to also have .bam.bai index files. Trying to run indexBam function on each"),"\n")
		indexBam(bamFiles)
		if(!all(file.exists(paste(bamFiles,".bai",sep="")))){
			stop("The bam files in UserDir are required to also have .bam.bai index files.")
		}else{
			if(verbose)cat(paste("Succesfully indexed all bamFiles in UserDir",UserDir),"\n")
		}
	}
	
	bamFilesList <- BamFileList(bamFiles)	#store all the .bam paths in a BamFile.
	
	#check that sequences in searchArea are actually found in the bam files
	header<-scanBamHeader(bamFiles)
	checkSeqNameExists<-function(bamHeader,requestedSeqNames){
		as.character(requestedSeqNames)%in%names(bamHeader[["targets"]])
	}
	if(!all(unlist(lapply(header,checkSeqNameExists,seqnames(searchArea))))){
		#not all searchArea requested seq-names found in bam files. Create nice error report and stop
		seqNotFoundErrors<-lapply(header,checkSeqNameExists,seqnames(searchArea))
		seqNotFounds<-vector()
		for(sampleName in names(seqNotFoundErrors)){
			seqNotFounds<-c(seqNotFounds,as.character(seqnames(searchArea)[!seqNotFoundErrors[[sampleName]]]))
		}
		stop(paste("The following seq name(s) not found in the bam files:",paste(sort(unique(seqNotFounds)),collapse=", ")))
	}
		
	#Loop through, open scanBam, store in GRList and then close each object in the BamFileList object.
	i<-1
	BamGRL<-GRangesList()
	for(bamName in names(bamFilesList)) {
		#Description
		bf <- bamFilesList[[bamName]]
		open(bf)
		if(verbose)cat(paste("Reading bam file",i,"with filename",basename(bamName)),"\n")	#Print information to the user
		bam<-scanBam(bf,param=param)
		#Description
		for(rangeName in names(bam)){
					
			ranges<-IRanges(
					start=bam[[rangeName]][["pos"]], #if NA values your in trouble. That means the read didnt map. Use better filter.
					width=cigarWidthAlongReferenceSpace(bam[[rangeName]][["cigar"]])	#bam[[rangeName]][["qwidth"]]	#bam[[rangeName]][["qwidth"]]
			)
			GRangeBam<-GRanges(
					seqnames = as.character(bam[[rangeName]][["rname"]]),   #Before "mrnm", now"rname"...
					ranges = ranges,
					strand = bam[[rangeName]][["strand"]],
					names=bam[[rangeName]][["qname"]],
					flag=bam[[rangeName]][["flag"]],
					cigar=bam[[rangeName]][["cigar"]],
					mapq=bam[[rangeName]][["mapq"]],
					mpos=bam[[rangeName]][["mpos"]],
					isize=bam[[rangeName]][["isize"]],
					seq=bam[[rangeName]][["seq"]],
					qual=bam[[rangeName]][["qual"]]	
			)
			
			#Store GRangeBam in BamGRL (which is the GRange List object)	
			if(basename(bamName)%in%names(BamGRL)){  #This way of merging the different chromosomes to the same GRangeObject is maybe not the best way. Later try to store them in a separate list, and then unlist before importing to GrangeBam
				BamGRL[[basename(bamName)]]<-c(BamGRL[[basename(bamName)]],GRangeBam)
			}
			else {
				BamGRL[[basename(bamName)]]<-GRangeBam
			}
		}
		if(verbose)cat(paste("stored",basename(bamName), "in BamGRL"),"\n")
		i<-1+i
		gc()
		close(bf)
	}
	return(BamGRL)
}


impBamGAL <- function(UserDir,searchArea,XStag=FALSE,verbose=TRUE){
	#Set parameters
	which <- searchArea #A GRanges, RangesList, RangedData, or missing object, from which a IRangesList instance will be constructed.
	what <- scanBamWhat() #A character vector naming the fields to return. scanBamWhat() returns a vector of available fields. Fields are described on the scanBam help page.
	flag <- scanBamFlag(isUnmappedQuery=FALSE)
	
	if(XStag){
		param<-ScanBamParam(flag=flag,which=which,what=what,tag="XS") #store ScanBamParam in param.
	}else{
		param<-ScanBamParam(flag=flag,which=which,what=what) #store ScanBamParam in param.		
	}
	#Point to correct directory and create a BamFileList object
	bamDir <- normalizePath(UserDir) #Point to the directory containing your Bam files and its respective bam.bai files.
	allFiles <- list.files(bamDir,full.names = TRUE) #list files in a folder.
	bamFiles <- allFiles[grep(".bam$",allFiles)] 	#list only the files ending in .bam .
	if(length(bamFiles) == 0)stop(paste("No bam files found in",bamDir))
	if(!all(file.exists(paste(bamFiles,".bai",sep="")))){
		if(verbose)cat(paste("The bam files in UserDir are required to also have .bam.bai index files. Trying to run indexBam function on each"),"\n")
		indexBam(bamFiles)
		if(!all(file.exists(paste(bamFiles,".bai",sep="")))){
			stop("The bam files in UserDir are required to also have .bam.bai index files.")
		}else{
			if(verbose)cat(paste("Succesfully indexed all bamFiles in UserDir",UserDir),"\n")
		}
	}
	bamFilesList <- BamFileList(bamFiles)	#store all the .bam paths in a BamFile.
	
	#check that sequences in searchArea are actually found in the bam files
	header<-scanBamHeader(bamFiles)
	checkSeqNameExists<-function(bamHeader,requestedSeqNames){
		as.character(requestedSeqNames)%in%names(bamHeader[["targets"]])
	}
	if(!all(unlist(lapply(header,checkSeqNameExists,seqnames(searchArea))))){
		#not all searchArea requested seq-names found in bam files. Create nice error report and stop
		seqNotFoundErrors<-lapply(header,checkSeqNameExists,seqnames(searchArea))
		seqNotFounds<-vector()
		for(sampleName in names(seqNotFoundErrors)){
			seqNotFounds<-c(seqNotFounds,as.character(seqnames(searchArea)[!seqNotFoundErrors[[sampleName]]]))
		}
		stop(paste("The following seq name(s) not found in the bam files:",paste(sort(unique(seqNotFounds)),collapse=", ")))
	}
	
	#Loop through, open scanBam, store in GRList and then close each object in the BamFileList object.
	BamGAL<-GAlignmentsList()
	i <- 1
	for(bamName in names(bamFilesList)) {
		#Description
		bf <- bamFilesList[[bamName]]
		open(bf)
		if(verbose)cat(paste("Reading bam file",i,"with filename",basename(bamName)),"\n")	#Print information to the user
		GappedAlign<-readGAlignmentsFromBam(bf,param=param)
		
		BamGAL[[basename(bamName)]]<-GappedAlign
		
		if(verbose)cat(paste("stored",basename(bamName), "in BamGAL"),"\n")
		gc()
		close(bf)
		i <- i +1
	}
	return(BamGAL)
}

impBcfGRL <- function(UserDir,searchArea,verbose=TRUE){
	
	#Set parameters
	which <- searchArea #A GRanges, RangesList, RangedData, or missing object, from which a IRangesList instance will be constructed.
	param <- ScanBcfParam(which=which)
	#Point to correct directory and create a BcfFileList object
	bcfDir <- normalizePath(UserDir) #Point to the directory containing your Bam files and its respective bam.bai files.
	allFiles <- list.files(bcfDir,full.names = TRUE) #list files in a folder.
	bcfFiles <- allFiles[grep(".bcf$",allFiles)] 	#list only the files ending in .bam .
	if(length(bcfFiles) == 0)stop(paste("No bcf files were found in",UserDir))
	
	#bcfFilesList <- BcfFileList(bcfFiles)	#store all the .bam paths in a BamFile.
	if(!all(file.exists(paste(bcfFiles,".bci",sep="")))){
		if(verbose)cat("Did not find bci files for all bcf files. Trying the indexBcf function obtain these","\n")
		for(bcfFile in bcfFiles){
			indexBcf(bcfFile)	
		}
		if(!all(file.exists(paste(bcfFiles,".bci",sep="")))){
			stop("The bcf files in UserDir are required to also have .bcf.bci index files. Run the indexBcf function in package Rsamtools on each bam file.")
		}
	}
	
	#Loop through, open scanBam, store in GRList and then close each object in the BamFileList object.
	BcfGRL<-GRangesList()
	for(i in 1:length(bcfFiles)) {

		bcf<-suppressWarnings(scanBcf(file=bcfFiles[i],param=param))
		
		#need to protect against empty bcf files
		if(length(bcf[["POS"]]) == 0){
			GRangeBcf<-GRanges(
					seqnames = vector(), 
					ranges = IRanges(
							start=vector(),
							width=vector() 
					),
					ref=vector(),					
					alt=vector(),
					qual=vector()
			)
			bcfName <- bcfFiles[i]
			BcfGRL[[basename(bcfName)]]<-GRangeBcf
			
		}else{ #if they are not empty we just proceed as usual
			
			ranges<-IRanges(
					start=bcf[["POS"]], #if NA values your in trouble. That means the read didnt map/Jesper
					width=1L #Width is set to "1" because its an Snp.
			)
			GRangeBcf<-GRanges(
					seqnames = as.character(bcf[["CHROM"]]), 
					ranges = ranges,
					ref=bcf[["REF"]],					
					alt=bcf[["ALT"]],
					qual=bcf[["QUAL"]]
			)
			#Store GRangeBam in BamGRL (which is the GRange List object)	
			bcfName <- bcfFiles[i]
			BcfGRL[[basename(bcfName)]]<-GRangeBcf
			if(verbose)cat(paste("stored",basename(bcfName), "in BcfGRL"),"\n")
			gc()
		}
	}
	return(BcfGRL)
}

#Function that wraps around BcfImpGRList, but returns a GRanges instead of GRangesList - since that seems to be more useful in many cases
impBcfGR <- function(UserDir, searchArea,verbose=TRUE){
	BcfGRList <- impBcfGRL(UserDir, searchArea, verbose)
	BcfGR<-do.call(c,unname(as.list(BcfGRList )))
	BcfGR<-unique(BcfGR)
	names(BcfGR)<-paste("chr",seqnames(BcfGR),"_",start(BcfGR),sep="")
	return(BcfGR)
}


realCigarPosition <- function(RleCigar,BpPos){

	#because of speed issues, checks are best performed outside this function.
	if(!class(RleCigar)=="Rle"){stop("class must be Rle")}

	e <- as.character(RleCigar)
	
	#changeVector
	v <- rep(0,length=length(e))
	names(v) <- e

	v[e=="M"] <- 1
	v[e=="I"] <- unlist(lapply(runLength(RleCigar)[runValue(RleCigar)=="I"],function(x){c(x+1,rep(1,x-1))}))
	#v[e=="D"] <- 0 #already zero 
	#v[e=="N"] <- 0 #already zero

	#sum all until interesting position
	cs <- cumsum(v)
	
	if(names(cs[BpPos])=="D"){retPos <- 0
	}else if(names(cs[BpPos])=="N"){retPos <- -1
	}else {
		retPos <- cs[BpPos]
		names(retPos) <- NULL
		if(retPos>sum(e=="M" | e=="I")){
			retPos <- -1 # the position went outside the read
		}
	}
	retPos
}



realCigarPositions <- function(RleCigar){
	#returns a vector that have order all positions to match with the cigar

	#because of speed issues, checks are best performed outside this function.
	if(!class(RleCigar)=="Rle"){stop("class must be Rle")}


	e <- as.character(RleCigar)	
	#make a new representation vector
	v <- rep(0,length=length(e))
	names(v) <- e

	v[e=="M"] <- 1
	v[e=="I"] <- unlist(lapply(runLength(RleCigar)[runValue(RleCigar)=="I"],function(x){c(x+1,rep(1,x-1))}))
	#v[e=="D"] <- 0
	#v[e=="N"] <- 0

	#sum all until interesting position
	cs <- cumsum(v)
	
	cs <- cs[!names(cs)=="D"]
	cs <- cs[!names(cs)=="N"]
	
	cs
}


realCigarPositionsList <- function(RleCigarList){

	#because of speed issues, checks are best performed outside this function.
	if(!class(RleCigarList)=="CompressedRleList"){stop("class must be Rle")}

	lapply(RleCigarList,
		function(RleCigar){
			e <- as.character(RleCigar)	
			#make a new representation vector
			v <- rep(0,length=length(e))
			names(v) <- e

			v[e=="M"] <- 1
			v[e=="I"] <- unlist(lapply(runLength(RleCigar)[runValue(RleCigar)=="I"],function(x){c(x+1,rep(1,x-1))}))
			#v[e=="D"] <- 0
			#v[e=="N"] <- 0

			#sum all until interesting position
			cs <- cumsum(v)
			
			cs <- cs[!names(cs)=="D"]
			cs <- cs[!names(cs)=="N"]
			
			cs
		}
	)
}




#create ldf object (i.e. for SnpAfList slot) based on a pre-defined list of Snps
getAlleleCount<-function(BamList, GRvariants,strand="nonStranded", verbose=TRUE){

	#if just one element of, make list (which is a convenient way of handling this input type)
	if(class(BamList)=="GAlignments"){BamList <- GAlignmentsList(BamList)}
	if(class(BamList)=="GRanges"){BamList <- GRangesList(BamList)}

	#check for strand name
	if(!class(strand)=="character"){stop("strand has to be of class character")}
	if(!length(strand)==1){stop("strand has to be of length 1")}
	if(!sum(strand %in% c("+","-","*","nonStranded"))>0){stop("strand parameter has to be either '+', '-', '*' or 'nonStranded' ")}

	#if the user sent in the GRangesList, take out only the unique entries.
	if(class(GRvariants)=="GRangesList"){	
	    GRvariants <- unique(unlist( GRvariants, use.names = FALSE)) #merge BcfGRL to one unique set of Snps
	}
	
	#checking that BamList format is ok (has to be done quite carefully for the list-class gappedAlignments 
	if(class(BamList)=="GRangesList"){
		warning("The use of GRangesList is not recommended. Use BamImpGAList or BamImpGAPList")
	}else if(class(BamList)=="GAlignmentsList"){ #assuming GAlignments or GAlignmentPairs and performing extra consistency checks
		#checks general for both GAlignments or GAlignmentPairs
		if(length(unique(unlist(lapply(BamList,class)))) != 1){
			stop(paste("BamList was given as a list containing entries of",length(unique(unlist(lapply(BamList,class)))),"different classes. This is not allowed"))
		}
		if(!unique(unlist(lapply(BamList,class)))%in%c("GAlignments","GAlignmentPairs")){
			stop(paste("BamList entries must be either of class	GAlignments or GAlignmentPairs not",unique(unlist(lapply(BamList,class)))))
		}
		
		#checks specific for GAlignments
		if(unique(unlist(lapply(BamList,class)))=="GAlignments"){
			if(!all(unlist(lapply(BamList,function(x){all(c("cigar","qwidth")%in%colnames(mcols(x)))})))){
				stop("BamList given as lists of GAlignments objects must contain mcols named qwidth and cigar")
			}
		}
		
		#checks specific for GAlignmentPairs
		if(unique(unlist(lapply(BamList,class)))=="GAlignmentPairs"){
			stop("BamList given as lists of GAlignmentPairs does not work yet (remove this when they do, and implement more consistency checks specific for these")
		}
	}else{ 
		stop("The class of the BamList has to be either GRangesList or a list with gappedAlignments or GAlignmentPairs")
	}
	
	#checking that GRvariants is ok
	if(class(GRvariants) != "GRanges")stop(paste("GRvariants must be of class GRanges, not",class(GRvariants)))
	if(length(GRvariants)== 0) stop("GRvariants was given as an empty GRanges object. There can be no Snps retrieved by getAlleleCount then")
	if(any(width(GRvariants)!=1))stop("GRvariants can contain only entries of width=1, corresponding to SNPs.")
	
	
	
	#checking that verbose is ok
	if(class(verbose) !="logical")stop(paste("verbose must be of class logical, not",class(verbose)))
	if(length(verbose) !=1)stop(paste("verbose must be of length 1, not",length(verbose)))
	
	#this added by /LF to protect against junction reads  (next 2 lines)
	if(class(BamList)=="GRangesList"){
		typicalReadLength<-as.numeric(names(table(width(ranges(unlist(BamList))))))[1] 
	} #this calculates the most common read-length 
	warnAgainstImpossibleReads<-FALSE
	
	#make row-names
	if(sum(grepl("chr",seqnames(GRvariants)))>0) { snpNames <- paste(seqnames(GRvariants),"_",start(GRvariants),sep="")
	}else{snpNames <- paste("chr",seqnames(GRvariants),"_",start(GRvariants),sep="")}
	
	ldf<-list()
	for(snp in snpNames){
		ldf[[snp]]<-matrix(rep(0,length(BamList)*5),nrow=length(BamList),ncol=5,dimnames=list(names(BamList),c("A","C","G","T","del")))
	}
	
	#find chromosomeLevels
	if(class(BamList)=="GRangesList"){
		chromosomeLevels <-levels(seqnames(unlist(BamList)))
	}else if(class(BamList)=="GAlignmentsList"){
		chromosomeLevels <- unique(unlist(lapply(BamList,function(x)levels(droplevels(runValue(seqnames(x)))))))
	}
	
	#use strand choice to only get reads from that strand
	if(!strand=="nonStranded"){BamList <- GAlignmentsList(mapply(function(x,y){x[y]},BamList,strand(BamList)==strand))}
	
	for(chr in chromosomeLevels){
		if(verbose)cat(chr,"\n")
		
		BamListchr <- GAlignmentsList(mapply(function(x,y){x[y]},BamList,seqnames(BamList)==chr))

		for(sample in names(BamListchr)){
			if(verbose)cat("sample ",sample,"\n")
			
			BamListHere<-BamListchr[[sample]]
			
			#this added by /LF to protect against junction reads (next 4 lines)
			if(class(BamList)=="GRangesList"){
				if(any((width(ranges(BamListHere)) > (typicalReadLength-3)) &	(width(ranges(BamListHere)) < (typicalReadLength+3)))){
					warnAgainstImpossibleReads<-TRUE
					BamListHere<-BamListHere[(width(ranges(BamListHere)) > (typicalReadLength-3)) &	(width(ranges(BamListHere)) < (typicalReadLength+3))]
				}	
			}
				
			#make interval tree over BamListHere
			subject <- BamListHere
			mcols(subject) <- NULL
			subject <- ranges(subject)				
			tree <- IntervalTree(subject)
						
			snpRange<-IRanges(start=start(GRvariants),width=1)
			
			if(!(length(snpRange)==0)){
				
				for(posNr in 1:length(snpRange)) {
					pos <- start(snpRange[posNr])
					if(verbose) cat(pos,"\t",posNr,"\n")					
					
					#find overlaps for this position
					hits <- findOverlaps(IRanges(start=pos,width=1),tree)
					TFrow <- countSubjectHits(hits) ==TRUE
					
					#get the reads here and the pos difference
					BamListThisPos<-BamListHere[TFrow]
					
					#ONLY for gappedAlignments 
					if(class(BamList)=="GAlignmentsList"){
						#Watch out! gaps can be present. use the cigar information 
						TFconsistentWidths <- width(BamListThisPos) == mcols(BamListThisPos)[["qwidth"]]
						TFcW <- TFconsistentWidths
						
						#premake boCount object for storing bpCount information	
						bpCount <- rep(NA,length(BamListThisPos))
						#count bp for "cigar all matched reads"
						bpCount[TFcW] <- pos - start(BamListThisPos[TFcW]) + 1
						
						if(sum(!TFcW)>0){

							bpTmp <- pos - start(BamListThisPos[!TFcW]) + 1
							
							#extract cigars
							cigar <- mcols(BamListThisPos[!TFcW])[["cigar"]]
							
							#make RleList for the cigars
							cl <- cigarToRleList(cigar)

							bpCount[!TFcW] <- mapply(realCigarPosition,cl,bpTmp)
							
							#Remove the reads that are not in the region of interest
							# -1 ,for skipped regions
							# 0 ,when position is a deletion
							TFr <- bpCount== -1 | bpCount== 0
							bpCount <- bpCount[!TFr]
							BamListThisPos <- BamListThisPos[!TFr]
							
							NrBpOnDeletion <- sum(bpCount== 0)

						}else{
							NrBpOnDeletion <- 0
						}
					
					}
						
					if(class(BamList)=="GRangesList"){
						if(length(BamListThisPos)>0){					
							bpCount <- start(GRvariants)[posNr] - start(ranges(BamListThisPos)) +1
							if(length(bpCount)>0){
								TFthisSnp <- (bpCount>0) & bpCount<width(ranges(BamListThisPos))
								bpCount<-bpCount[TFthisSnp]
								BamListThisPos<-BamListThisPos[TFthisSnp]
								if(any(!TFthisSnp)){ 
									warnAgainstImpossibleReads<-TRUE
								}
							}
						}	
					}
					
					if(length(bpCount)>0){
						readNucleotide<-as.character(subseq(mcols(BamListThisPos)[["seq"]], start=bpCount, width=1))
						ldf[[snpNames[posNr]]][sample,"A"]<-sum(readNucleotide%in%"A")
						ldf[[snpNames[posNr]]][sample,"C"]<-sum(readNucleotide%in%"C")
						ldf[[snpNames[posNr]]][sample,"G"]<-sum(readNucleotide%in%"G")
						ldf[[snpNames[posNr]]][sample,"T"]<-sum(readNucleotide%in%"T")
						ldf[[snpNames[posNr]]][sample,"del"] <- NrBpOnDeletion
					}
				}
			}
		}
	}
	#this added by /LF to protect against junction reads (next 1 line)
	if(warnAgainstImpossibleReads)warning("The BamList contained reads with indicated widths longer than the length of the actual contained sequence")
	
	ldf
}




#investigate BamList sequence data for possible heterozygotes. Returns a vector of positions, which can be used with extractFromUserDefinedSnps to retrieve counts 
scanForHeterozygotes<-function(BamList, minimumReadsAtPos = 20, maximumMajorAlleleFrequency = 0.9, minimumBiAllelicFrequency = 0.9,maxReads=15000, verbose=TRUE){

	#if just one element of, make list (which is a convenient way of handling this input type)
	if(class(BamList)=="GAlignments"){BamList <- GAlignmentsList(BamList)}
	if(class(BamList)=="GRanges"){BamList <- GRangesList(BamList)}

	if(class(minimumReadsAtPos) !="numeric")stop(paste("minimumReadsAtPos must be of class numeric, not",class(minimumReadsAtPos)))
	if(length(minimumReadsAtPos) !=1)stop(paste("minimumReadsAtPos must be of length 1, not",length(minimumReadsAtPos)))
	if(class(maximumMajorAlleleFrequency) !="numeric")stop(paste("maximumMajorAlleleFrequency must be of class numeric, not",class(maximumMajorAlleleFrequency)))
	if(length(maximumMajorAlleleFrequency) !=1)stop(paste("maximumMajorAlleleFrequency must be of length 1, not",length(maximumMajorAlleleFrequency)))
	if(maximumMajorAlleleFrequency < 0 | maximumMajorAlleleFrequency >1)stop("maximumMajorAlleleFrequency must be between 0 and 1") 
	if(class(minimumBiAllelicFrequency) !="numeric")stop(paste("minimumBiAllelicFrequency must be of class numeric, not",class(minimumBiAllelicFrequency)))
	if(length(minimumBiAllelicFrequency) !=1)stop(paste("minimumBiAllelicFrequency must be of length 1, not",length(minimumBiAllelicFrequency)))
	if(minimumBiAllelicFrequency < 0 | maximumMajorAlleleFrequency >1)stop("minimumBiAllelicFrequency must be between 0 and 1") 
	if(class(verbose) !="logical")stop(paste("verbose must be of class logical, not",class(verbose)))
	if(length(verbose) !=1)stop(paste("verbose must be of length 1, not",length(verbose)))
	
	#checking that BamList format is ok (has to be done quite carefully for the list-class gappedAlignments 
	if(class(BamList)=="GRangesList")stop("The use of GRangesList is not recommended. Use BamImpGAList or BamImpGAPList")
	if(class(BamList)!="GAlignmentsList")stop("The class of the BamList has to be a GAlignmentsList")

	#checks specific for GAlignments
	if(unique(unlist(lapply(BamList,class)))=="GAlignments"){
		if(!all(unlist(lapply(BamList,function(x){all(c("cigar","qwidth")%in%colnames(mcols(x)))})))){
			stop("BamList given as GAlignmentsLists of GAlignments objects must contain mcols named qwidth and cigar")
		}
	}

	#check that we dont create a too large and memory-consuming matrix
	if(sum(unlist(lapply(BamList,length)) > maxReads)>0){stop("you may consume too much memory. If there is plenty of memory, then increase maxReads to allow more reads")}
	
	RangedData <- GRanges()
	chromosomeLevels <- unique(unlist(lapply(BamList,function(x){levels(droplevels(runValue(seqnames(x))))})))
	
	for(chr in chromosomeLevels){
		if(verbose)cat(paste("Investigating chromosome",chr),"\n")

		BamListchr <- GAlignmentsList(mapply(function(x,y){x[y]},BamList,seqnames(BamList)==chr))


		for(sample in 1:length(BamListchr)){
			if(verbose)cat(paste("Investigating sample",sample),"out of", length(BamListchr),"\n")
			
			#extract samples
			BamListHere<-BamListchr[[sample]]
			
			if(!(length(BamListHere)==0)){

				#then iterate over all reads 
				cigarRle<-cigarToRleList(mcols(BamListHere)[,"cigar"])
				toKeep <- realCigarPositionsList(cigarRle)
				
				seq<-mcols(BamListHere)[,"seq"]
				charList <- strsplit(as.character(seq),"")
				
				start <- start(BamListHere)-min(start(BamListHere))+1

				nw <- max(end(BamListHere))-min(start(BamListHere)) +2

				#populate matrix
				new <- matrix(NA, nrow=max(end(BamListHere))-min(start(BamListHere))+1, ncol=length(BamListHere))
				for(i in 1:ncol(new)){
					new[start[i]:(start[i]+length(toKeep[[i]])-1),i] <- charList[[i]][toKeep[[i]]]
				
				}

				#set rownames
				rownames(new)<-as.character(1:(nw-1))

				new<-new[apply(!is.na(new),1,sum) > minimumReadsAtPos,]
				
				if(!nrow(new)==0) {	

					#tabulate countsPerPosition (cpp)
					cpp<-apply(new,1,table)
				
					TFl <-unlist(lapply(cpp,
						function(x){
							if(length(x)>1){
								
								MajorAlleleFrequency<-x[order(x,decreasing=TRUE)[1]]  / sum(x)
								if(MajorAlleleFrequency < maximumMajorAlleleFrequency){
									MinorAlleleFrequency<-x[order(x,decreasing=TRUE)[2]]  / sum(x)
									if((MinorAlleleFrequency + MajorAlleleFrequency) > minimumBiAllelicFrequency){
										TRUE
									}else{FALSE}
								}else{FALSE}		
							}else{FALSE}
						}
					))
					if(!all(!TFl)) {		
						GR <- GRanges(ranges=IRanges(start=(as.numeric(names(cpp[TFl]))+min(start(BamListHere))-1),width=1),seqnames=chr)
						RangedData <- c(RangedData,GR)
					}
				}
			}
		}
	}
	#merge from all individuals
	RangedData<-unique(RangedData)
	
	#Add a Snp name based on position
	if(!(length(RangedData)==0)) {
		names(RangedData)<-paste("chr",seqnames(RangedData),"_",start(RangedData),sep="")	
		
	}
	return(RangedData)
}






decorateWithGenes <- function(x,genesInRegion,xlim,ylim,chromosome){
	
	#check the input variables
	if(class(xlim) != "integer")xlim<-as.numeric(xlim)
	if(class(xlim) != "numeric")stop(paste("xlim must be of class numeric, not",class(xlim)))
	if(length(xlim) != 2)stop(paste("xlim must be of length 2, not",length(xlim)))
	if(class(ylim) != "integer")ylim<-as.numeric(ylim)
	if(class(ylim) != "numeric")stop(paste("ylim must be of class numeric, not",class(ylim)))
	if(length(ylim) != 2)stop(paste("ylim must be of length 2, not",length(ylim)))
	if(class(chromosome) != "character")stop(paste("chromosome must be of class character, not",class(chromosome)))
	if(length(chromosome) != 1)stop(paste("chromosome must be of length 1, not",length(chromosome)))
	if(!chromosome%in%unique(seqnames(genesInRegion))){
		if(sub("^chr","",chromosome)%in%unique(seqnames(genesInRegion))){
			chromosome<-sub("^chr","",chromosome)				
		}else{
			stop(paste("chromosome",chromosome,"was not found amongst the seqnames of the genesInRegion object:",paste(sort(unique(seqnames(genesInRegion))),collapse=", ")))
		}
	}
	
	
	#only work with genes on current chromosome
	genesInRegion<-genesInRegion[seqnames(genesInRegion)==chromosome]
	
	#calculate how many "rows" have to be made available
	maxCoverage <- max(coverage(genesInRegion)[[chromosome]])
	
	#loop over all unique genes, drawing them as specified
	uniqueGenes<-unique(mcols(genesInRegion)[["Symbol"]])
	
	for(i in 1:length(uniqueGenes)){
		#getting the name of the gene and the height on the Y-axis 
		genesymbol <- uniqueGenes[i]
		yPos <- ylim[1] + (i-1)%%maxCoverage * ((ylim[2] - ylim[1]) / maxCoverage)
		
		#this block checks for double instances and just arbitrarily take the first (typically miRNAs with two locations)
		if(sum(mcols(genesInRegion)[["Symbol"]]%in%genesymbol) > 1){
			geneData<-genesInRegion[which(mcols(genesInRegion)[["Symbol"]]%in%genesymbol)[1],]
		}else{
			geneData<-genesInRegion[which(mcols(genesInRegion)[["Symbol"]]%in%genesymbol),]	
		}
		
		#draw and label
		start <- max(c(xlim[1]-(xlim[2]-xlim[1])/10,start(geneData)))
		end <- min(c(xlim[2]+(xlim[2]-xlim[1])/10,end(geneData)))
		lines(x = c(start,end), y=c(yPos,yPos),lwd=2)
		text(x = start + (end-start)/2 , y = yPos+(ylim[2]-ylim[1])/6, label=genesymbol,cex=0.8)
	}
}





decorateWithExons <- function(x,exonsInRegion,xlim,ylim,chromosome){
	
	#check the input variables
	if(class(exonsInRegion)!="GRanges")stop(paste("exonsInRegion must be of class GRanges, not",class(exonsInRegion)))
	if(!"tx_name"%in%colnames(mcols(exonsInRegion)))stop("exonsInRegion must contain an mcol variable named 'tx_name'")
	if(class(xlim) != "integer")xlim<-as.numeric(xlim)
	if(class(xlim) != "numeric")stop(paste("xlim must be of class numeric, not",class(xlim)))
	if(length(xlim) != 2)stop(paste("xlim must be of length 2, not",length(xlim)))
	if(class(ylim) != "integer")ylim<-as.numeric(ylim)
	if(class(ylim) != "numeric")stop(paste("ylim must be of class numeric, not",class(ylim)))
	if(length(ylim) != 2)stop(paste("ylim must be of length 2, not",length(ylim)))
	if(class(chromosome) != "character")stop(paste("chromosome must be of class character, not",class(chromosome)))
	if(length(chromosome) != 1)stop(paste("chromosome must be of length 1, not",length(chromosome)))
	if(!chromosome%in%unique(seqnames(exonsInRegion))){
		if(sub("^chr","",chromosome)%in%unique(seqnames(exonsInRegion))){
			chromosome<-sub("^chr","",chromosome)				
		}else{
			stop(paste("chromosome",chromosome,"was not found amongst the seqnames of the exonInRegion object:",paste(sort(unique(seqnames(exonsInRegion))),collapse=", ")))
		}
	}
	#only work with exons on current chromosome
	exonsInRegion<-exonsInRegion[seqnames(exonsInRegion)==chromosome]
	
	
	#calculate how many "rows" have to be made available (corresponding to the number of unique transcripts in exonsInRegion
	uniqueGenes<-unique(unlist(as.list(mcols(exonsInRegion)[["tx_name"]])))
	maxCoverage <- length(uniqueGenes)
	

	for(i in 1:length(uniqueGenes)){
		#getting the name of the gene and the height on the Y-axis 
		tx_name <- uniqueGenes[i]
		yPos <- ylim[1] + (i-1)%%maxCoverage * ((ylim[2] - ylim[1]) / maxCoverage)
		
		#extracting and drawing all exons for each transcript
		exonsInTranscript<-which(sapply(as.list(mcols(exonsInRegion)[["tx_name"]]),function(x){tx_name%in%x}))
		for(exonInTranscript in exonsInTranscript){
			x1<-start(exonsInRegion[exonInTranscript])
			x2<-end(exonsInRegion[exonInTranscript])
			lines(x=c(x1,x2),y=c(yPos,yPos),lwd=2)
		}
		
		
		#draw a thin connecting line for each transcript	
		end<-max(end(exonsInRegion[exonsInTranscript])) 
		start<-min(start(exonsInRegion[exonsInTranscript]))
		lines(x=c(start,end),y=c(yPos,yPos),lwd=0.5)		
		
		#label with the tx_name
		start<-max(c(start,xlim[1]))	
		end<-min(c(end,xlim[2]))
		text(x = start + (end-start)/2 , y = yPos+(ylim[2]-ylim[1])/6, label=tx_name,cex=0.8)
		
	}
}

getGenesFromAnnotation <- function(OrgDb, GR, leftFlank=1000,rightFlank=1000, getUCSC=FALSE, verbose=FALSE) {
	#checks
	if(class(OrgDb)!="OrgDb")stop(paste("OrgDb must of class OrgDb, not",class(OrgDb)))
	
	if(class(GR)!="GRanges")stop(paste("GR must of class GRanges, not",class(GR)))
	
	if(!class(leftFlank)%in%c("numeric"))stop(paste("leftFlank must be of class numeric, not:",class(leftFlank)))
	if(length(leftFlank)!=1)stop(paste("leftFlank must be of length 1, not:",length(leftFlank)))
	if(leftFlank<0)stop(paste("leftFlank must be equal to or larger than 0"))
	
	if(!class(rightFlank)%in%c("numeric"))stop(paste("rightFlank must be of class numeric, not:",class(rightFlank)))
	if(length(rightFlank)!=1)stop(paste("rightFlank must be of length 1, not:",length(rightFlank)))
	if(rightFlank<0)stop(paste("rightFlank must be equal to or larger than 0"))

	if(!class(getUCSC)%in%c("logical"))stop(paste("getUCSC must be of class logical, not:",class(getUCSC)))
	if(length(getUCSC)!=1)stop(paste("getUCSC must be of length 1, not:",length(getUCSC)))
	
	if(!"UCSCKG"%in%columns(OrgDb)){
		if(verbose)message("Unable to retrieve UCSCKG column from OrgDb. Omitting")
		getUCSC<-FALSE
	}
	
	if(!class(verbose)%in%c("logical"))stop(paste("verbose must be of class logical, not:",class(verbose)))
	if(length(verbose)!=1)stop(paste("verbose must be of length 1, not:",length(verbose)))
	
	#remove chr in seqnames
	seqLevels <- sub("^chr", "", seqlevels(GR))
	seqlevels(GR) <- seqLevels
	
	#pre-filtering to local region +/- 1MB (for speed purposes) 
	startFilter<- max(c(1,start(range(GR)) - 10^6))
	endFilter<- end(range(GR)) + 10^6
	colsFilter<-c("CHR","CHRLOC","CHRLOCEND","SYMBOL")
	sFilter<-suppressWarnings(select(OrgDb, keys=seqLevels, cols=colsFilter , keytype="CHR"))
	symbolsToGet<- sFilter[abs(sFilter[,"CHRLOC"]) > startFilter & abs(sFilter[,"CHRLOCEND"]) < endFilter & !is.na(sFilter[,"CHRLOCEND"]) & !is.na(sFilter[,"CHRLOC"]),"SYMBOL"] 
	
	
	#then create an annGR for genes in range
	if(getUCSC){
		cols <- c("SYMBOL","CHR","CHRLOC","CHRLOCEND","ENSEMBL","UCSCKG")
	}else{
		cols <- c("SYMBOL","CHR","CHRLOC","CHRLOCEND","ENSEMBL")
	}
	s <- suppressWarnings(select(OrgDb, keys=symbolsToGet, cols=cols, keytype="SYMBOL"))
	
	
	#remove Symbols with NAs
	TFminusStrand <- s[["CHRLOC"]] <0
	TFplusStrand <- s[["CHRLOC"]] >0
	sNoNas <- s[c(which(TFminusStrand),which(TFplusStrand)),]
	
	#make Strand vector
	TFminusStrand2 <- sNoNas[["CHRLOC"]] <0
	strand <- rep("+",length=(dim(sNoNas)[1] ))
	
	strand[TFminusStrand2] <- "-"
	
	#make start and end vector 
	sNonNegative <- sNoNas
	sNonNegative[TFminusStrand2,c("CHRLOC","CHRLOCEND")] <- -sNonNegative[TFminusStrand2,c("CHRLOC","CHRLOCEND")]
	start <- sNonNegative[["CHRLOC"]]
	end <- sNonNegative[["CHRLOCEND"]]
	
	#make seqnames
	seqnames <- sNonNegative[["CHR"]]
	
	#make the annGR containing all genes
	if(getUCSC){	
		annGR <- GRanges(
				seqnames=Rle(seqnames),
				ranges=IRanges(start,end),
				strand= Rle(strand),
				Symbol=sNonNegative[["SYMBOL"]],
				Ensembl=sNonNegative[["ENSEMBL"]],
				UCSCKG=sNonNegative[["UCSCKG"]]
		)
	}else{
		annGR <- GRanges(
				seqnames=Rle(seqnames),
				ranges=IRanges(start,end),
				strand= Rle(strand),
				Symbol=sNonNegative[["SYMBOL"]],
				Ensembl=sNonNegative[["ENSEMBL"]]
		)
		
	}
	
	#check that all levels in GR exist in annGR, if not exclude these levels
	if(sum(!levels(seqnames(GR)) %in% levels(seqnames(annGR)))>0){
		TFkeepLevels <- levels(seqnames(GR)) %in% levels(seqnames(annGR))
		seqlevels(GR, force=FALSE) <- seqlevels(GR)[TFkeepLevels]
	
	}

	#check that all levels in annGR exist in GR, if not exclude these levels
	if(sum(!levels(seqnames(annGR)) %in% levels(seqnames(GR)))>0){
		TFkeepLevels <- levels(seqnames(annGR)) %in% levels(seqnames(GR))
		seqlevels(annGR, force=TRUE) <- seqlevels(annGR)[TFkeepLevels]	
	}

	#the seqlevels comes in different orders. This will give the correct order.
	seqlevels(GR) <- seqlevels(annGR)

	#find overlaps between annGR and Snps incl. flank region
	GenesInRegion <- subsetByOverlaps(annGR,GR) #put in flankSize here when you have time ;)
	seqlengths(GenesInRegion) <- seqlengths(GR)
	GenesInRegion
}




getGenesVector <- function(OrgDb, GR,verbose=FALSE){
	if(verbose){cat("start gene extraction\n")}
	GenesInRegion <- getGenesFromAnnotation(OrgDb, GR, leftFlank=1000,rightFlank=1000,verbose=FALSE)	

	seqlevels(GR) <- seqlevels(GenesInRegion)

	#remove duplicate symbol names
	#if same symbol merge regions
	symbolList <- unique(mcols(GenesInRegion)[["Symbol"]])
	newGenesInRegion <- GRanges()
	
	#check if GenesInRegion is zero
	if(!length(GenesInRegion)==0){
		for(i in 1:length(symbolList)){
			symbol <- symbolList[i]
			TF <- mcols(GenesInRegion)[["Symbol"]]==symbol

			G <- GRanges(
				seqnames=unique(seqnames(GenesInRegion[TF])),
				ranges=IRanges(min(start(GenesInRegion[TF])),max(end(GenesInRegion[TF]))),
				strand=unique(strand(GenesInRegion[TF])))
			mcols(G) <- unique(mcols(GenesInRegion[TF])[,"Symbol",drop=FALSE])

		newGenesInRegion <- c(newGenesInRegion,G) 
		}
	}

	#half-vectorized solution
	h <- findOverlaps(newGenesInRegion,GR)
	symbolVec <- vector() 
	for(i in 1:length(GR)){
		symbolVec[i] <- paste(mcols(newGenesInRegion[queryHits(h)[subjectHits(h)==i]])[["Symbol"]], collapse=",")
	}
	#set NAs where appropriate
	symbolVec[symbolVec==""] <- NA
	#return list with symbols
	symbolVec
	
	#return list with symbols
	symbolVec
}

getExonsFromAnnotation <- function(TxDb, GR,leftFlank=1000,rightFlank=1000,verbose=FALSE) {
	
	#checks
	if(class(TxDb)!="TranscriptDb")stop(paste("GR must of class TranscriptDb, not",class(TxDb)))
	
	if(class(GR)!="GRanges")stop(paste("GR must of class GRanges, not",class(GR)))
		
	if(!class(leftFlank)%in%c("numeric"))stop(paste("leftFlank must be of class numeric, not:",class(leftFlank)))
	if(length(leftFlank)!=1)stop(paste("leftFlank must be of length 1, not:",length(leftFlank)))
	if(leftFlank<0)stop(paste("leftFlank must be equal to or larger than 0"))
	
	if(!class(rightFlank)%in%c("numeric"))stop(paste("rightFlank must be of class numeric, not:",class(rightFlank)))
	if(length(rightFlank)!=1)stop(paste("rightFlank must be of length 1, not:",length(rightFlank)))
	if(rightFlank<0)stop(paste("rightFlank must be equal to or larger than 0"))
	
	if(!class(verbose)%in%c("logical"))stop(paste("verbose must be of class logical, not:",class(verbose)))
	if(length(verbose)!=1)stop(paste("verbose must be of length 1, not:",length(verbose)))
	
	#remove chr in seqnames for GR
	seqLevels <- sub("^chr", "", seqlevels(GR))
	seqlevels(GR) <- seqLevels

	seqlevels(TxDb, force=TRUE) <- paste("chr",names(seqlengths(GR)),sep="")

	#Get all exons from the active chromosomes
	#By creating a GRanges from TxDb
	annGR <- exons(TxDb, columns=c("exon_id","tx_name"))

	#remove chr in seqnames for annGR
	seqLevels <- sub("^chr", "", seqlevels(annGR))
	seqlevels(annGR) <- seqLevels
	
	#check that all levels in GR exist in annGR, if not exclude these levels
	if(sum(!levels(seqnames(GR)) %in% levels(seqnames(annGR)))>0){
		TFkeepLevels <- levels(seqnames(GR)) %in% levels(seqnames(annGR))
		seqlevels(GR, force=FALSE) <- seqlevels(GR)[TFkeepLevels]
		
	}

	#check that all levels in annGR exist in GR, if not exclude these levels
	if(sum(!levels(seqnames(annGR)) %in% levels(seqnames(GR)))>0){
		TFkeepLevels <- levels(seqnames(annGR)) %in% levels(seqnames(GR))
		seqlevels(annGR, force=TRUE) <- seqlevels(annGR)[TFkeepLevels]	
	}

	#the seqlevels comes in different orders. This will give the correct order.
	seqlevels(GR) <- seqlevels(annGR)

	#add flanking regions
	lf <- flank(GR,leftFlank,start=TRUE)
	rf <- flank(GR,rightFlank,start=FALSE)
	GR <- c(lf,GR,rf)
	start<-max(c(1,min(start(GR))))	
	end<-max(end(GR))	
	
	#speed-increasing coarse subset to +/- 1MB before extracting (because smaller annGR is faster and we have to access several times)
	annGR<-annGR[start(ranges(annGR))> max(c(1,start-10^6)) & end(ranges(annGR)) < (end+10^6)]
	
	#extract names of all transcripts with exons in plotting window
	tx_names<-unique(unlist(as.list(mcols(annGR[start(ranges(annGR))>start & end(ranges(annGR)) < end])[["tx_name"]])))
	
	#retrieve all transcript annotation, including potential off-plot tails and heads
	ExonsInRegion<-annGR[sapply(mcols(annGR)[["tx_name"]],function(x){any(tx_names%in%x)})]
	
	seqlengths(ExonsInRegion) <- seqlengths(GR)
	ExonsInRegion

}


getExonsVector <- function(TxDb, GR,verbose=FALSE){
	if(verbose){cat("start exon extraction\n")}

	ExonsInRegion <- getExonsFromAnnotation(TxDb, GR,leftFlank=0,rightFlank=0,verbose=verbose)

	seqlevels(GR) <- seqlevels(ExonsInRegion)

	#half-vectorized solution
	h <- findOverlaps(ExonsInRegion,GR)
	ExonVec <- vector() 
	for(i in 1:length(GR)){
		ExonVec[i] <- paste(mcols(ExonsInRegion[queryHits(h)[subjectHits(h)==i]])[["exon_id"]], collapse=",")
	}
	#set NAs where appropriate
	ExonVec[ExonVec==""] <- NA
	#return list with symbols
	ExonVec

}


getTranscriptsFromAnnotation <- function(TxDb, GR,leftFlank=1000,rightFlank=1000,verbose=FALSE) {
	
	#checks
	if(class(TxDb)!="TranscriptDb")stop(paste("GR must of class TranscriptDb, not",class(TxDb)))
	
	if(class(GR)!="GRanges")stop(paste("GR must of class GRanges, not",class(GR)))
	
	if(!class(leftFlank)%in%c("numeric"))stop(paste("leftFlank must be of class numeric, not:",class(leftFlank)))
	if(length(leftFlank)!=1)stop(paste("leftFlank must be of length 1, not:",length(leftFlank)))
	if(leftFlank<0)stop(paste("leftFlank must be equal to or larger than 0"))
	
	if(!class(rightFlank)%in%c("numeric"))stop(paste("rightFlank must be of class numeric, not:",class(rightFlank)))
	if(length(rightFlank)!=1)stop(paste("rightFlank must be of length 1, not:",length(rightFlank)))
	if(rightFlank<0)stop(paste("rightFlank must be equal to or larger than 0"))
	
	if(!class(verbose)%in%c("logical"))stop(paste("verbose must be of class logical, not:",class(verbose)))
	if(length(verbose)!=1)stop(paste("verbose must be of length 1, not:",length(verbose)))
	
	
	
	#remove chr in seqnames for GR
	seqLevels <- sub("^chr", "", seqlevels(GR))
	seqlevels(GR) <- seqLevels

	seqlevels(TxDb, force=TRUE) <- paste("chr",names(seqlengths(GR)),sep="")


	#Get all exons from the active chromosomes
	#By creating a GRanges from TxDb
	annGR <- transcripts(TxDb)


	#remove chr in seqnames for annGR
	seqLevels <- sub("^chr", "", seqlevels(annGR))
	seqlevels(annGR) <- seqLevels
	
	#check that all levels in GR exist in annGR, if not exclude these levels
	if(sum(!levels(seqnames(GR)) %in% levels(seqnames(annGR)))>0){
		TFkeepLevels <- levels(seqnames(GR)) %in% levels(seqnames(annGR))
		seqlevels(GR, force=FALSE) <- seqlevels(GR)[TFkeepLevels]
		
	}

	#check that all levels in annGR exist in GR, if not exclude these levels
	if(sum(!levels(seqnames(annGR)) %in% levels(seqnames(GR)))>0){
		TFkeepLevels <- levels(seqnames(annGR)) %in% levels(seqnames(GR))
		seqlevels(annGR, force=TRUE) <- seqlevels(annGR)[TFkeepLevels]	
	}

	#the seqlevels comes in different orders. This will give the correct order.
	seqlevels(GR) <- seqlevels(annGR)

	#add flanking regions
	lf <- flank(GR,leftFlank,start=TRUE)
	rf <- flank(GR,rightFlank,start=FALSE)
	GR <- c(lf,GR,rf)
	#find overlaps between annGR and Snps incl. flank region
	TxInRegion <- subsetByOverlaps(annGR,GR) #put in flankSize here when you have time ;)
	seqlengths(TxInRegion) <- seqlengths(GR)
	TxInRegion
}

getTranscriptsVector <- function(TxDb, GR,verbose=FALSE){
	if(verbose){cat("start transcript extraction\n")}

	TxInRegion <- getTranscriptsFromAnnotation(TxDb, GR,leftFlank=0,rightFlank=0)

	seqlevels(GR) <- seqlevels(TxInRegion)

	#half-vectorized solution
	h <- findOverlaps(TxInRegion,GR)
	TxVec <- vector() 
	for(i in 1:length(GR)){
		TxVec[i] <- paste(mcols(TxInRegion[queryHits(h)[subjectHits(h)==i]])[["tx_id"]], collapse=",")
	}
	#set NAs where appropriate
	TxVec[TxVec==""] <- NA
	#return list with symbols
	TxVec

}

getCDSFromAnnotation <- function(TxDb, GR,leftFlank=1000,rightFlank=1000,verbose=FALSE) {
	#CDS are the coding regions that do not only code for proteins, but other also other types like RNA.

	#checks
	if(class(TxDb)!="TranscriptDb")stop(paste("GR must of class TranscriptDb, not",class(TxDb)))
	
	if(class(GR)!="GRanges")stop(paste("GR must of class GRanges, not",class(GR)))
	
	if(!class(leftFlank)%in%c("numeric"))stop(paste("leftFlank must be of class numeric, not:",class(leftFlank)))
	if(length(leftFlank)!=1)stop(paste("leftFlank must be of length 1, not:",length(leftFlank)))
	if(leftFlank<0)stop(paste("leftFlank must be equal to or larger than 0"))
	
	if(!class(rightFlank)%in%c("numeric"))stop(paste("rightFlank must be of class numeric, not:",class(rightFlank)))
	if(length(rightFlank)!=1)stop(paste("rightFlank must be of length 1, not:",length(rightFlank)))
	if(rightFlank<0)stop(paste("rightFlank must be equal to or larger than 0"))
	
	if(!class(verbose)%in%c("logical"))stop(paste("verbose must be of class logical, not:",class(verbose)))
	if(length(verbose)!=1)stop(paste("verbose must be of length 1, not:",length(verbose)))
	
	#remove chr in seqnames for GR
	seqLevels <- sub("^chr", "", seqlevels(GR))
	seqlevels(GR) <- seqLevels

	seqlevels(TxDb, force=TRUE) <- paste("chr",names(seqlengths(GR)),sep="")


	#Get all exons from the active chromosomes
	#By creating a GRanges from TxDb
	annGR <- cds(TxDb)

	
	#remove chr in seqnames for annGR
	seqLevels <- sub("^chr", "", seqlevels(annGR))
	seqlevels(annGR) <- seqLevels
	
	#check that all levels in GR exist in annGR, if not exclude these levels
	if(sum(!levels(seqnames(GR)) %in% levels(seqnames(annGR)))>0){
		TFkeepLevels <- levels(seqnames(GR)) %in% levels(seqnames(annGR))
		seqlevels(GR, force=FALSE) <- seqlevels(GR)[TFkeepLevels]
		
	}

	#check that all levels in annGR exist in GR, if not exclude these levels
	if(sum(!levels(seqnames(annGR)) %in% levels(seqnames(GR)))>0){
		TFkeepLevels <- levels(seqnames(annGR)) %in% levels(seqnames(GR))
		seqlevels(annGR, force=TRUE) <- seqlevels(annGR)[TFkeepLevels]	
	}

	#the seqlevels comes in different orders. This will give the correct order.
	seqlevels(GR) <- seqlevels(annGR)

	#add flanking regions
	lf <- flank(GR,leftFlank,start=TRUE)
	rf <- flank(GR,rightFlank,start=FALSE)
	GR <- c(lf,GR,rf)
	#find overlaps between annGR and Snps incl. flank region
	CDSInRegion <- subsetByOverlaps(annGR,GR) #put in flankSize here when you have time ;)
	seqlengths(CDSInRegion) <- seqlengths(GR)
	CDSInRegion
	
}



getCDSVector <- function(TxDb, GR,verbose=FALSE){
	if(verbose){cat("start CDS extraction\n")}

	CDSInRegion <- getCDSFromAnnotation(TxDb, GR,leftFlank=0,rightFlank=0)

	seqlevels(GR) <- seqlevels(CDSInRegion)

	#half-vectorized solution
	h <- findOverlaps(CDSInRegion,GR)
	CDSVec <- vector() 
	for(i in 1:length(GR)){
		CDSVec[i] <- paste(mcols(CDSInRegion[queryHits(h)[subjectHits(h)==i]])[["cds_id"]], collapse=",")
	}
	#set NAs where appropriate
	CDSVec[CDSVec==""] <- NA
	#return list with symbols
	CDSVec
}

getAnnotationDataFrame <- function(GR,strand="+",annotationType=NULL,OrgDb=NULL, TxDb=NULL,verbose=FALSE)
{
	#main checks
	if(sum(!(annotationType %in% c("gene","exon","cds","transcript")))>0){
		stop("annotationType must be one or more of these arguments 'gene','exon','cds','transcript'")}

	if(is.null(OrgDb) & is.null(TxDb)){
		stop("at least one of parameters OrgDb or TxDb must be used")}
	
	#nr of columns for return df
	ncol <- 0
	if(!is.null(OrgDb)){ncol <- ncol + 1}
	if(!is.null(TxDb)){ncol <- ncol + (length(annotationType) -1)}	 
	
	#return dataframe
	df <- data.frame(row.names=1:length(GR))
	
	#set strand
	strand(GR) <- strand
	
	#extract annotation
	if(!is.null(OrgDb)){
		if("gene"%in%annotationType){
			gene <- getGenesVector(OrgDb=OrgDb, GR, verbose=verbose)
			df[["symbol"]] <- gene
		}
		if(is.null(annotationType)){
			gene <- getGenesVector(OrgDb=OrgDb, GR, verbose=verbose)
			df[["symbol"]] <- gene
		}
	}
	if(!is.null(TxDb)){
		if("exon"%in%annotationType){
			df[["exon_id"]] <-getExonsVector(TxDb=TxDb,GR, verbose=verbose)
		}
		if("transcript"%in%annotationType){
			df[["tx_id"]] <- getTranscriptsVector(TxDb=TxDb,GR, verbose=verbose)
		}
		if("cds"%in%annotationType){
			df[["cds_id"]] <- getCDSVector(TxDb=TxDb,GR, verbose=verbose)
		}

		if(is.null(annotationType)) {
			df[["exon_id"]] <-getExonsVector(TxDb=TxDb,GR, verbose=verbose)
			df[["tx_id"]] <- getTranscriptsVector(TxDb=TxDb,GR, verbose=verbose)
			df[["cds_id"]] <- getCDSVector(TxDb=TxDb,GR, verbose=verbose)
		}
	}
	df
}

getDefaultMapBiasExpMean <- function(alleleCountList){

	l <- lapply(alleleCountList,function(x){
			ap <- apply(x,2,sum)
			char <- names(sort(ap,decreasing=TRUE))[1:2]
			
			v <- rep(0,length(colnames(x)))
			v[colnames(x) %in% char] <- 0.5
			v
		}
	)

	MapBiasExpMean <- matrix(unlist(l),byrow=TRUE,nrow=length(alleleCountList),ncol=5,dimnames=list(c(names(alleleCountList)),colnames(alleleCountList[[1]]))) # alleleCountList[[1]] assumes that in each list the colnames are the same.
	MapBiasExpMean
}

getDefaultMapBiasExpMean3D <- function(alleleCountList){

	MapBiasExpMean <- getDefaultMapBiasExpMean(alleleCountList)
	#make 3D array	
	MapBiasExpMean3D <- array(NA,c(length(alleleCountList),length(unlist(unique(lapply(alleleCountList,rownames)))),5)) #empty array
	for(i in 1:length(unlist(unique(lapply(alleleCountList,rownames))))) {
		MapBiasExpMean3D[,i,] <- MapBiasExpMean
	}
	MapBiasExpMean3D
}



extractReferenceAllele <- function(GR,path){
	#GR contains ranges of positions to be retrieved
	#path contains the path to fasta reference used for alignment i.e. ("hg19.fa")
	
	#check if there exist an index
	#Check to be written

	#remove chromosome name
	if(!(grepl("chr",seqlevels(GR)))) {
		seqLevels <- paste("chr", seqlevels(GR),sep="")
		seqlevels(GR) <- seqLevels
	}

	fl <- FaFile(path)

	#open file
	open(fl)
	#scan fasta
	res <- scanFa(fl,GR)
	#clsoe file
	close(fl)
	res

}

getAreaFromGeneNames <- function(genesymbols, OrgDb, leftFlank=1000,rightFlank=1000, verbose=TRUE){
	
	#start up sets
	if(!class(genesymbols)%in%c("character"))stop(paste("genesymbols must be of class character, not:",class(genesymbols)))
	
	if(!class(OrgDb)%in%c("OrgDb"))stop(paste("OrgDb must be of class OrgDb, not:",class(OrgDb)))
	
	if(!class(leftFlank)%in%c("numeric"))stop(paste("leftFlank must be of class numeric, not:",class(leftFlank)))
	if(length(leftFlank)!=1)stop(paste("leftFlank must be of length 1, not:",length(leftFlank)))
	if(leftFlank<0)stop(paste("leftFlank must be equal to or larger than 0"))
	
	if(!class(rightFlank)%in%c("numeric"))stop(paste("rightFlank must be of class numeric, not:",class(rightFlank)))
	if(length(rightFlank)!=1)stop(paste("rightFlank must be of length 1, not:",length(rightFlank)))
	if(rightFlank<0)stop(paste("rightFlank must be equal to or larger than 0"))
	
	if(!class(verbose)%in%c("logical"))stop(paste("verbose must be of class logical, not:",class(verbose)))
	if(length(verbose)!=1)stop(paste("verbose must be of length 1, not:",length(verbose)))
	
	#retrieving data
	colsFilter<-c("CHR","CHRLOC","CHRLOCEND","SYMBOL")
	s <-suppressWarnings(select(OrgDb, keys=genesymbols, cols=colsFilter , keytype="SYMBOL"))
	
	missing<-genesymbols[!genesymbols%in%s[,"SYMBOL"]]
	if(verbose & length(missing)>0){
		cat(paste("Did not find information on these",length(missing),"genes:",paste(missing,collapse=", ")),"\n") 	
	}else{
		if(verbose)cat("Found all requested genes in annotation","\n")
	}
	
	searchArea<-GRanges(
			seqnames = paste("chr",s[,"CHR"],sep=""),
			ranges = IRanges(
					abs(s[,"CHRLOC"])-leftFlank,
					abs(s[,"CHRLOCEND"])+rightFlank
			)
	)
	
	searchArea<-reduce(searchArea)
	
	return(searchArea)
}







getSnpIdFromLocation <- function(GR, SNPloc, verbose=TRUE) {
	
	if(class(GR)!="GRanges")stop(paste("GR must of class GRanges, not",class(GR)))
	if(class(SNPloc)!="SNPlocs")stop(paste("SNPlocs must of class SNPlocs, not",class(SNPloc)))
	if(!exists("getSNPlocs"))stop("There must exist a function called getSNPlocs, available from the SNPlocs.Hsapiens.dbSNP.xxxxxxx package. Try to load such package.")
	
	# add chr to seqnames if not present
	if(length(grep("^chr",seqnames(GR)))!=length(GR)){
		if(length(grep("^chr",seqnames(GR)))!=0)stop("seqnames must all begin with 'chr'. In the GR it seemed that some, but not all, seqnames began with chr. Please correct manually")
		seqnames<-paste("chr",seqnames(GR),sep="")
		seqlevels(GR)<-as.character(unique(seqnames))
		seqnames(GR)<-seqnames
	}
	
	#changing chr to ch to adapt to SNPloc
	if(length(grep("^chr",seqnames(GR)))==length(GR)){
		seqnames<-seqnames(GR)
		seqnames<-sub("^chr","ch",seqnames(GR))
		seqlevels(GR)<-as.character(unique(seqnames))
		seqnames(GR)<-seqnames
	}
	
	
	SNPlocThisChr<-getSNPlocs(seqlevels(GR), as.GRanges=TRUE, caching=FALSE)
	
	overlaps<-findOverlaps(GR, SNPlocThisChr) 
	
	if(verbose)cat(paste("Replacing position-based SNP name with rs-ID for",length(overlaps),"SNP(s)"),"\n")
	
	#replace name in GR
	for(i in 1:length(overlaps)){
		snp<-paste("rs",mcols(SNPlocThisChr[subjectHits(overlaps[i])])[,"RefSNP_id"],sep="")
		names(GR)[queryHits(overlaps[i])] <-snp
	}

	#change back to chr from ch
	seqlevels(GR) <- sub("^ch","chr",seqlevels(GR))

	return(GR)
}


