#'@include ASEset-class.R
NULL

#' Import Bam
#' 
#' Imports a specified genomic region from a bam file using a GenomicRanges
#' object as search area.
#' 
#' These functions are wrappers to import bam files into R and store them into
#' either GRanges, GAlignments or GappedAlignmentpairs objects.
#' 
#' It is recommended to use the impBamGAL() which takes information of gaps
#' into account. It is also possible to use the other variants as well, but
#' then pre-filtering becomes important because gapped, intron-spanning reads
#' will cause problems. This is because the GRanges objects can not handle if
#' gaps are present and will then give a wrong result when calculating the
#' allele (SNP) count table.
#' 
#' If the sequence data is strand-specific you may want to set XStag=TRUE. The
#' strand specific information will then be stored in the meta columns with
#' column name 'XS'.
#' 
#' @name import-bam
#' @rdname import-bam
#' @aliases import-bam impBamGAL impBamGRL
#' @param UserDir The relative or full path of folder containing bam files.
#' @param searchArea A \code{GenomicRanges object} that contains the regions of
#' interest
#' @param XStag Setting \code{XStag=TRUE} stores the strand specific
#' information in the mcols slot 'XS'
#' @param verbose Setting \code{verbose=TRUE} gives details of procedure during
#' function run.
#' @return \code{impBamGRL} returns a GRangesList object containing the RNA-seq
#' reads in the region defined by the \code{searchArea} argument.
#' \code{impBamGAL} returns a list with GAlignments objects containing the
#' RNA-seq reads in the region defined by the \code{searchArea} argument.
#' \code{funImpBamGAPL} returns a list with GappedAlignmentPairs object
#' containing the RNA-seq reads in the region defined by the \code{searchArea}
#' argument.
#' @note A typical next step after the import of bam files is to obtain SNP
#' information. This can be done either with the \code{\link{impBcfGRL}} and
#' \code{\link{getAlleleCounts}} functions. Alternatively the
#' \code{\link{scanForHeterozygotes}} function provides R-based functionality
#' for identifying heterozygote coding SNPs.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{impBcfGRL}} for importing Bcf
#' files.  }
#' @keywords bam import
#' @examples
#' 
#' #Declare searchArea
#' searchArea <- GRanges(seqnames=c("17"), ranges=IRanges(79478301,79478361))
#' 
#' #Relative or full path  
#' pathToFiles <- system.file("extdata/ERP000101_subset", package="AllelicImbalance")
#' 
#' reads <- impBamGAL(pathToFiles,searchArea,verbose=FALSE)
#' 
#' @export impBamGRL
#' @export impBamGAL
NULL


#' @rdname import-bam
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

#' @rdname import-bam
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
	BamGAL<- list()
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
	BamGAL <- GAlignmentsList(BamGAL)

	return(BamGAL)
}

#' Import Bcf Selection
#' 
#' Imports a selection of a bcf file or files specified by a GenomicRanges
#' object as search area.
#' 
#' A wrapper to import bcf files into R in the form of GenomicRanges objects.
#' 
#' @name import-bcf
#' @rdname import-bcf
#' @aliases import-bcf impBcfGRL impBcfGR
#' @param UserDir The relative or full path of folder containing bam files.
#' @param searchArea A \code{GenomicRanges} object that contains the regions of
#' interest
#' @param verbose Setting \code{verbose=TRUE} gives details of the procedure
#' during function run.
#' @return \code{BcfImpGRList} returns a GRangesList object.  \code{BcfImpGR}
#' returns one GRanges object of all unique entries from one or more bcf files.
#' @note Make sure there is a complementary index file \code{*.bcf.bci} for
#' each bcf file in \code{UserDir}. If there is not, then the functions
#' \code{impBcfGRL} and \code{impBcfGR} will try to create them.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{impBamGRL}} for importing bam files
#' \item The \code{\link{getAlleleCounts}} for how to get allele(SNP) counts
#' \item The \code{\link{scanForHeterozygotes}} for how to find possible
#' heterozygote positions }
#' @keywords bcf import
#' @examples
#' 
#' #Declare searchArea
#' searchArea <- GRanges(seqnames=c("17"), ranges=IRanges(79478301,79478361))
#' 
#' #Relative or full path  
#' pathToFiles <- system.file("extdata/ERP000101_subset", package="AllelicImbalance")
#' 
#' #import
#' reads <- impBcfGRL(pathToFiles, searchArea, verbose=FALSE)
#' 
#' 
#' @export impBcfGRL
#' @export impBcfGR
NULL

#' @rdname import-bcf
impBcfGRL <- function(UserDir,searchArea=NULL,verbose=TRUE){
	
	#Set parameters
	if(is.null(searchArea)){ 
		param <- ScanBcfParam()
	}else{ 
		param <- ScanBcfParam(which=searchArea)
	}
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

#' @rdname import-bcf
impBcfGR <- function(UserDir, searchArea=NULL,verbose=TRUE){
	BcfGRList <- impBcfGRL(UserDir, searchArea, verbose)
	BcfGR<-do.call(c,unname(as.list(BcfGRList )))
	BcfGR<-unique(BcfGR)
	names(BcfGR)<-paste("chr",seqnames(BcfGR),"_",start(BcfGR),sep="")
	return(BcfGR)
}


#' realCigarPosition
#' 
#' From a GAlignments calculate the real corresponding position for each read
#' based on its cigar.
#' 
#' The main intention for these functions are to be the internal functions for
#' \code{scanForHeterozygotes} and \code{getAlleleCount}.
#' 
#' @name cigar-utilities
#' @rdname cigar-utilities
#' @aliases cigar-utilities realCigarPosition realCigarPositions
#' realCigarPositionsList
#' @param RleCigar An \code{Rle} containing cigar information
#' @param RleCigarList An \code{RleList} containing cigar information
#' @param BpPos the absolute position on the chromosome of interest
#' @return \code{realCigarPosition} returns the new position
#' \code{realCigarPositions} returns a vector with the corrected positions to
#' be subsetted from a read.  \code{realCigarPositionsList} returns a list
#' where each element i a vector with the corrected positions to be subsetted
#' from a read.
#' @author Jesper R. Gadin
#' @seealso \itemize{ \item The \code{\link{scanForHeterozygotes}} which is a
#' function to find possible heterozygote sites in a
#' \code{\link[GenomicAlignments]{GAlignmentsList}} object }
#' @keywords internal
#' @examples
#' 
#'   RleCigarList <-  cigarToRleList("3M4I93M")
#'   BpPos <- 5
#' 
#'   newPos <- realCigarPosition(RleCigar=RleCigarList[[1]], BpPos)
#'   newPositions <- realCigarPositions(RleCigar=RleCigarList[[1]])
#'   newPositionsList <- realCigarPositionsList(RleCigarList=RleCigarList)
#' @export realCigarPosition
#' @export realCigarPositions
#' @export realCigarPositionsList
NULL

#' @rdname cigar-utilities
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



#' @rdname cigar-utilities
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


#' @rdname cigar-utilities
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

#' snp count data
#' 
#' Given the positions of known SNPs, this function returns allele counts from
#' a BamGRL object
#' 
#' This function is used to retrieve the allele counts from specified positions
#' in a set of RNA-seq reads. The \code{BamList} argument will typically have
#' been created using the \code{impBamGAL} function on bam-files. The
#' \code{GRvariants} is either a GRanges with user-specified locations or else
#' it is generated through scanning the same bam-files as in \code{BamList} for
#' heterozygote locations (e.g. using \code{scanForHeterozygotes}). The
#' GRvariants will currently only accept locations having width=1,
#' corresponding to SNPs.  In the \code{strand} argument, specifying
#' 'nonStranded' is the same as retrieving the sum count of '+' and '-' reads
#' (and '*' unknown in case these are found in the bam file). 'nonStranded' is
#' the default behaviour and can be used when the RNA-seq experiments strand
#' information is not available.
#' 
#' @param BamList A \code{GAlignmentsList object} or \code{GRangesList object}
#' containing data imported from a bam file
#' @param GRvariants A \code{GRanges object} that contains positions of SNPs to
#' retrieve
#' @param strand A length 1 \code{character} with value 'nonStranded', '+',
#' '-', or '*'.  This argument determines if \code{getAlleleCounts} will
#' retrieve counts from all reads, or only from reads marked as '+', '-' or '*'
#' (unknown), respectively.
#' @param return.type "list" or "array"
#' @param verbose Setting \code{verbose=TRUE} makes function more talkative
#' @return \code{getAlleleCounts} returns a list of several data.frame objects,
#' each storing the count data for one SNP.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{scanForHeterozygotes}} which is a
#' function to find possible heterozygote sites in a
#' \code{\link[GenomicAlignments]{GAlignmentsList}} object }
#' @keywords SNP count
#' @examples
#' 
#' 	#load example data
#' 	data(reads)
#' 	data(GRvariants)
#' 
#' 	#set seqlevels in reads equal to seqlevels in GRvariants
#' 	seqlevels(reads) <- "17"
#' 
#' 	#get counts at the three positions specified in GRvariants
#' 	alleleCount <- getAlleleCounts(BamList=reads,GRvariants,
#' 		strand="nonStranded")
#' 	
#' 	#if the reads had contained stranded data, these two calls would 
#' 	#have given the correct input objects for getAlleleCounts
#' 	alleleCountPlus <- getAlleleCounts(BamList=reads,GRvariants,
#' 		strand="+")
#' 	alleleCountMinus <- getAlleleCounts(BamList=reads,GRvariants,
#' 		strand="-")
#' 	
#' 
#' @export getAlleleCounts
getAlleleCounts<-function(BamList, GRvariants, strand="nonStranded", return.type="list", verbose=TRUE){

	if(!class(BamList)%in%c("GAlignments","GAlignmentsList")){stop("BamList has to be of class GAlignments or GAlignmnetsList\n")}
	#if just one element of, make list (which is a convenient way of handling this input type)
	if(class(BamList)=="GAlignments"){BamList <- GAlignmentsList(BamList)}

	#check for strand name
	if(!class(strand)=="character"){stop("strand has to be of class character")}
	if(!length(strand)==1){stop("strand has to be of length 1")}
	if(!sum(strand %in% c("+","-","*","nonStranded"))>0){stop("strand parameter has to be either '+', '-', '*' or 'nonStranded' ")}

	#if the user sent in the GRangesList for GRvariants, take out only the unique entries.
	if(class(GRvariants)=="GRangesList"){	
	    GRvariants <- unique(unlist( GRvariants, use.names = FALSE)) #merge BcfGRL to one unique set of Snps
	}

	#check that seqlevels are the same	
	if(!identical(seqlevels(BamList), seqlevels(GRvariants))){stop("!identical(seqlevels(BamList), seqlevels(GRvariants))\n")}

	
	#checking that GRvariants is ok
	if(class(GRvariants) != "GRanges")stop(paste("GRvariants must be of class GRanges, not",class(GRvariants)))
	if(length(GRvariants)== 0) stop("GRvariants was given as an empty GRanges object. There can be no Snps retrieved by getAlleleCount then")
	if(any(width(GRvariants)!=1))stop("GRvariants can contain only entries of width=1, corresponding to SNPs.")
	
	
	
	#checking that verbose is ok
	if(class(verbose) !="logical")stop(paste("verbose must be of class logical, not",class(verbose)))
	if(length(verbose) !=1)stop(paste("verbose must be of length 1, not",length(verbose)))
	
	
	#make row-names
	if(sum(grepl("chr",seqnames(GRvariants)))>0) { snpNames <- paste(seqnames(GRvariants),"_",start(GRvariants),sep="")
	}else{snpNames <- paste("chr",seqnames(GRvariants),"_",start(GRvariants),sep="")}
	
	dimnames= list(snpNames,names(BamList),c("A","C","G","T"))
	ar1 <- array(NA,c(length(GRvariants),length(BamList),4),dimnames=dimnames) #empty array that handles only four nucleotides + one del columns
	
	#use strand choice to only get reads from that strand
	if(!strand=="nonStranded"){BamList <- GAlignmentsList(mapply(function(x,y){x[y]},BamList,strand(BamList)==strand))}

	for(j in 1:length(names(BamList))){
		sample <- names(BamList)[j]	
		if(verbose)cat("sample ",sample,"\n")

		gal <-BamList[[j]]
	
		nuclpiles <- pileLettersAt(mcols(gal)[,"seq"], seqnames(gal), start(gal), cigar(gal), GRvariants)

		#fill array
		nstr <- strsplit(as.character(nuclpiles),"")
		for(k in 1:length(GRvariants)){
			ar1[k,j,] <- c(sum(nstr[[k]]%in%"A"),sum(nstr[[k]]%in%"C"),sum(nstr[[k]]%in%"G"),sum(nstr[[k]]%in%"T")) #del will always be 0. Could have set it to NA, but then it makes problem further  down in the chain of functions...			
		}

		
	}

	#check return.type argument	
	if(return.type=="list"){
		alleleCountList <- list()
		for(i in 1:nrow(ar1)){
			mat <- ar1[i,,]
			if(class(mat)=="numeric"){
					mat <- t(mat)
					colnames(mat) <- dimnames[[3]]
				}else{
					colnames(mat) <- dimnames[[3]]
				}
			rownames(mat) <- dimnames[[2]] 
			alleleCountList[[i]] <- mat
		}
		names(alleleCountList) <- dimnames[[1]]	
		alleleCountList
	}else if(return.type=="array"){
		ar1
	}else{cat("return.type unknown\n Nothing will be returned from function!")}
}





#create ldf object (i.e. for SnpAfList slot) based on a pre-defined list of Snps
#old version


#' snp count data
#' 
#' Given the positions of known SNPs, this function returns allele counts from
#' a BamGRL object
#' 
#' This function is used to retrieve the allele counts from specified positions
#' in a set of RNA-seq reads. The \code{BamList} argument will typically have
#' been created using the \code{impBamGAL} function on bam-files. The
#' \code{GRvariants} is either a GRanges with user-specified locations or else
#' it is generated through scanning the same bam-files as in \code{BamList} for
#' heterozygote locations (e.g. using \code{scanForHeterozygotes}). The
#' GRvariants will currently only accept locations having width=1,
#' corresponding to SNPs.  In the \code{strand} argument, specifying
#' 'nonStranded' is the same as retrieving the sum count of '+' and '-' reads
#' (and '*' unknown in case these are found in the bam file). 'nonStranded' is
#' the default behaviour and can be used when the RNA-seq experiments strand
#' information is not available.
#' 
#' @param BamList A \code{GAlignmentsList object} or \code{GRangesList object}
#' containing data imported from a bam file
#' @param GRvariants A \code{GRanges object} that contains positions of SNPs to
#' retrieve
#' @param strand A length 1 \code{character} with value 'nonStranded', '+',
#' '-', or '*'.  This argument determines if \code{getAlleleCounts2} will
#' retrieve counts from all reads, or only from reads marked as '+', '-' or '*'
#' (unknown), respectively.
#' @param verbose Setting \code{verbose=TRUE} makes function more talkative
#' @return \code{getAlleleCounts2} returns a list of several data.frame
#' objects, each storing the count data for one SNP.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{scanForHeterozygotes}} which is a
#' function to find possible heterozygote sites in a
#' \code{\link[GenomicAlignments]{GAlignmentsList}} object }
#' @keywords SNP count
#' @examples
#' 
#' 	#load example data
#' 	data(reads)
#' 	data(GRvariants)
#' 
#' 	#set seqlevels in reads equal to seqlevels in GRvariants
#' 	seqlevels(reads) <- "17"
#' 
#' 	#get counts at the three positions specified in GRvariants
#' 	alleleCount <- getAlleleCounts2(BamList=reads,GRvariants,
#' 		strand="nonStranded")
#' 	
#' 	#if the reads had contained stranded data, these two calls would 
#' 	#have given the correct input objects for getAlleleCounts2
#' 	alleleCountPlus <- getAlleleCounts2(BamList=reads,GRvariants,
#' 		strand="+")
#' 	alleleCountMinus <- getAlleleCounts2(BamList=reads,GRvariants,
#' 		strand="-")
#' 	
#' 
#' @export getAlleleCounts2
getAlleleCounts2<-function(BamList, GRvariants,strand="nonStranded", verbose=TRUE){

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


#' scanForHeterozygotes
#' 
#' Identifies the positions of SNPs found in BamGR reads.
#' 
#' This function scans all reads stored in a \code{GAlignmentsList} for
#' possible heterozygote positions. The user can balance the sensitivity of the
#' search by modifying the minimumReadsAtPos, maximumMajorAlleleFrequency and
#' minimumBiAllelicFrequency arguments.
#' 
#' @param BamList A \code{GAlignmentsList object}
#' @param minimumReadsAtPos minimum number of reads required to call a SNP at a
#' given position
#' @param maximumMajorAlleleFrequency maximum frequency allowed for the most
#' common allele. Setting this parameter lower will minimise the SNP calls
#' resulting from technical read errors, at the cost of missing loci with
#' potential strong ASE
#' @param minimumBiAllelicFrequency minimum frequency allowed for the first and
#' second most common allele. Setting a Lower value for this parameter will
#' minimise the identification of loci with three or more alleles in one
#' sample. This is useful if sequencing errors are suspected to be common.
#' @param maxReads max number of reads of one list-element allowed
#' @param verbose logical indicating if process information should be displayed
#' @return \code{scanForHeterozygotes} returns a GRanges object with the SNPs
#' for the BamList object that was used as input.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{getAlleleCounts}} which is a
#' function that count the number of reads overlapping a site.  }
#' @keywords scan SNP heterozygote
#' @examples
#' 
#' 	data(reads)
#' 	s <- scanForHeterozygotes(reads,verbose=FALSE)
#' 
#' @export scanForHeterozygotes
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






#' Map Bias
#' 
#' an allele frequency list
#' 
#' This function will assume there is no bias that comes from the mapping of
#' reads, and therefore create a matrix with expected frequency of 0.5 for each
#' allele.
#' 
#' @aliases getDefaultMapBiasExpMean getDefaultMapBiasExpMean3D
#' @param alleleCountList A \code{GRangesList object} containing read
#' information
#' @return \code{getDefaultMapBiasExpMean} returns a matrix with a default
#' expected mean of 0.5 for every element.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords mapping bias
#' @examples
#' 
#' 	#load example data
#' 	data(ASEset)
#' 	#access SnpAfList
#' 	alleleCountList <- alleleCounts(ASEset)
#' 	#get default map bias exp mean
#' 	matExpMean <- getDefaultMapBiasExpMean(alleleCountList)
#' 
#' 
#' @export getDefaultMapBiasExpMean
getDefaultMapBiasExpMean <- function(alleleCountList){

	l <- lapply(alleleCountList,function(x){
			ap <- apply(x,2,sum)
			char <- names(sort(ap,decreasing=TRUE))[1:2]
			
			v <- rep(0,length(colnames(x)))
			v[colnames(x) %in% char] <- 0.5
			v
		}
	)

	MapBiasExpMean <- matrix(unlist(l),byrow=TRUE,nrow=length(alleleCountList),ncol=4,dimnames=list(c(names(alleleCountList)),colnames(alleleCountList[[1]]))) # alleleCountList[[1]] assumes that in each list the colnames are the same.
	MapBiasExpMean
}

getDefaultMapBiasExpMean3D <- function(alleleCountList){

	MapBiasExpMean <- getDefaultMapBiasExpMean(alleleCountList)
	#make 3D array	
	MapBiasExpMean3D <- array(NA,c(length(alleleCountList),length(unlist(unique(lapply(alleleCountList,rownames)))),4)) #empty array
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
#' 	#load example data
#' 	data(ASEset)
#' 	
#' 	#get counts at the three positions specified in GRvariants
#' 	library(org.Hs.eg.db )
#' 	searchArea<-getAreaFromGeneNames(c("PAX8","TLR7"), org.Hs.eg.db)	
#' 
#' @export getAreaFromGeneNames
getAreaFromGeneNames <- function(genesymbols, OrgDb, leftFlank=0,rightFlank=0,na.rm=FALSE, verbose=TRUE){
	
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
	
	#remove NAs
	if(na.rm==TRUE){
		s <- s[!(is.na(s[,"CHR"]) | is.na(s[,"CHRLOC"]) | is.na(s[,"CHRLOCEND"])),]
	}else{
		warningGenes <- s[is.na(s[,"CHR"]) | is.na(s[,"CHRLOC"]) | is.na(s[,"CHRLOCEND"]),]
		if(!nrow(warningGenes)==0){
			warning(paste(warningGenes[,"SYMBOL"],"had NAs","\n"),"you better remove these genes from your 'genesymbols'")
		}
	}

	strand <- rep("+",nrow(s))
	TFstrand <- s[,"CHRLOC"] < 0
	strand[TFstrand] <- "-"

	searchArea<-GRanges(
			seqnames = paste("chr",s[,"CHR"],sep=""),
			ranges = IRanges(
					abs(s[,"CHRLOC"])-leftFlank,
					abs(s[,"CHRLOCEND"])+rightFlank
			),
			strand=strand,
			symbol = s[["SYMBOL"]]
	)
	
	searchArea<-reduce(searchArea, with.revmap=TRUE)
	l <- lapply(searchArea$revmap,function(x){
			paste(unique(s[x,"SYMBOL"]),collapse=",")		
		}
	)
	mcols(searchArea)[,"symbol"] <- unlist(l)
	#remove the mapping column
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
#' 	#load example data
#' 	data(ASEset)
#' 	
#' 	#get counts at the three positions specified in GRvariants
#' 	if(require(SNPlocs.Hsapiens.dbSNP.20120608)){
#'		updatedGRanges<-getSnpIdFromLocation(rowData(ASEset), SNPlocs.Hsapiens.dbSNP.20120608)
#'		rowData(ASEset)<-updatedGRanges
#'	}
#' 
#' @export getSnpIdFromLocation
getSnpIdFromLocation <- function(GR, SNPloc, return.vector=FALSE, verbose=TRUE) {
	
	if(class(GR)!="GRanges")stop(paste("GR must of class GRanges, not",class(GR)))
	if(class(SNPloc)!="SNPlocs")stop(paste("SNPlocs must of class SNPlocs, not",class(SNPloc)))
	if(!exists("getSNPlocs"))stop("There must exist a function called getSNPlocs, available from the SNPlocs.Hsapiens.dbSNP.xxxxxxx package. Try to load such package.")
	
	# add chr to seqnames if not present
	if(length(grep("^chr",seqnames(GR)))!=length(GR)){
		if(length(grep("^chr",seqnames(GR)))!=0)stop("seqnames must all begin with 'chr'. In the GR it seemed that some, but not all, seqnames began with chr. Please correct manually")
		seqlevels(GR)<-paste("chr",seqlevels(GR),sep="")
		#seqnames(GR)<-seqnames
	}
	
	#changing chr to ch to adapt to SNPloc
	if(length(grep("^chr",seqnames(GR)))==length(GR)){
		#seqnames<-seqnames(GR)
		seqlevels(GR)<-sub("^chr","ch",seqlevels(GR))
		#seqlevels(GR)<-as.character(unique(seqnames))
		#seqlevels(GR) <- levels(seqnames)
		#seqnames(GR)<-seqnames
	}
	
	
	SNPlocThisChr<-getSNPlocs(seqlevels(GR), as.GRanges=TRUE, caching=FALSE)

	seqlevels(SNPlocThisChr,force=TRUE) <- seqlevels(GR) 		
#	seqlengths(GR) <- seqlengths(SNPlocThisChr)

	overlaps<-findOverlaps(GR, SNPlocThisChr) 
	
	if(verbose)cat(paste("Replacing position-based SNP name with rs-ID for",length(overlaps),"SNP(s)"),"\n")
	
	#replace name in GR
	#for(i in 1:length(overlaps)){
	#	snp<-paste("rs",mcols(SNPlocThisChr[subjectHits(overlaps[i])])[,"RefSNP_id"],sep="")
	#	names(GR)[queryHits(overlaps[i])] <-snp
	#}

	snp<-paste("rs",mcols(SNPlocThisChr[subjectHits(overlaps)])[,"RefSNP_id"],sep="")
	names(GR)[queryHits(overlaps)] <-snp
		
	#change back to chr from ch
	seqlevels(GR) <- sub("^ch","chr",seqlevels(GR))
	
	if(return.vector){names(GR)
	}else{return(GR)}

}

#Deprecated functions
getAlleleCount <- function()
{
	    .Deprecated("getAlleleCounts")
    ## use new function, or remainder of myOldFunc
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
#' @param afraction the fractions in matrix for each snp
#' @param arank the ranks for each snp
#' @param acounts the counts for each snp
#' @param amainVec, vector for each plots main
#' @param ... for simpler generics when extending function
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{ASEset}} class which the barplot
#' function can be called up on.  }
#' @keywords barplot
#' @examples
#' 
#' 	a <- ASEset
#' 	strand <- "+"
#' 	acounts <-  alleleCounts(a,strand=strand)
#' 	arank <-  arank(a,strand=strand)
#' 	afraction <- fraction(a, strand=strand)
#' 	amainVec <- rep("",nrow(a))
#' 
#' 	name <- rownames(a)[1]
#' 
#' 	barplotLatticeFraction(identifier=name, afraction, arank, amainVec) 
#' 	barplotLatticeCounts(identifier=name,  acounts, arank, amainVec) 
#' 
#' @export barplotLatticeFraction
#' @export barplotLatticeCounts
#' 
NULL

#' @rdname barplot-lattice-support
barplotLatticeFraction <- function(identifier, afraction, arank, amainVec, ... ){
#afraction 
#arank
	gargs <- list(...)

	a.r <- arank[[identifier]][1:2]	
	a.f <- afraction[,identifier]
	a.f2 <- 1-a.f
	values <- vector()
	for (i in 1:length(a.f)){values <- c(values,a.f[i],a.f2[i])}
	allele <- rep(a.r,length(a.f))

	#sample <- names(a.f)
	sample <- vector()
	for (i in 1:length(a.f)){sample <- c(sample,rownames(afraction)[i],rownames(afraction)[i])}
	df <- data.frame(values=values,sample=sample,allele=allele)

	TFna <- is.na(df$values)
	df$values[TFna] <- 0 # 0.5 + 0.5 -> 1
	na <-rep("no",length(values)) 
	na[TFna] <- "yes"
	df <- cbind(df,na)

	df$sample <- factor(df$sample,levels=unique(df$sample))
	#df$sample <- factor(df$sample,levels=rownames(df))


	#Grah params
	my_cols <- c("green", "red")

	#set default values 
	parset<- list()

	scales = list(rot=c(90,0))
	gargs$deAnnoPlot <- FALSE

	if(gargs$deAnnoPlot){

		parset<- list(
			 layout.widths = list(
			 left.padding = 0,
			 axis.left = 0,
			 ylab.axis.padding =0,
			 right.padding = 0,
			 axis.right = 0
		))

		scales = list(y=list(at=NULL,labels=NULL),rot=c(90,0))

		#tmp
		gargs$ylab <- ""
		gargs$xlab <- ""
	}

 	 b <- barchart(values~sample,
		 #horiz=FALSE,
		 group=allele,
		 data=df,
		 col = my_cols,
		 origin=0,
		 #auto.key=list(points = FALSE, rectangles = TRUE,space="top",size=2,cex=0.8),
		 stack=TRUE,
		 scales = scales,
		 main=amainVec,
		 ylab=gargs$ylab,
		 xlab=gargs$xlab,
		 par.settings=parset
		 #box.ratio=2,
		 #abbreviate=TRUE
	)

	b

}

#' @rdname barplot-lattice-support
barplotLatticeCounts <- function(identifier, acounts, arank, amainVec, ...){
	
	gargs <- list(...)

	#implodeList(gargs)

	a.m <- amainVec[identifier]
	a.r <- arank[[identifier]][1:2]	
	a.c <- acounts[[identifier]][,a.r,drop=FALSE]

	values <- as.vector(t(a.c))
	allele <- rep(colnames(a.c),nrow(a.c))

	sample <- vector()
	for (i in 1:nrow(a.c)){sample <- c(sample,rownames(a.c)[i],rownames(a.c)[i])}
	df <- data.frame(values=values,sample=sample,allele=allele)

	#to get right order in barchart
	df$sample <- factor(df$sample,levels=unique(df$sample))
	#df$values[is.na(df$values)] <- 0 #doesnt work

	###
	#Grah params
	###
	#set default values 
	parset<- list()
	
	scales = list(rot=c(90,0))
	gargs$deAnnoPlot <- FALSE

	#potentially override default settings with trellis settings
	if(gargs$deAnnoPlot){

		parset<- list(
			 layout.widths = list(
			 left.padding = 0,
			 axis.left = 0,
			 ylab.axis.padding =0,
			 right.padding = 0,
			 axis.right = 0
		))


		scales = list(y=list(at=NULL,labels=NULL),rot=c(90,0))
	}

	b <- barchart(values~sample,
	 horiz=FALSE,
	 origin=0,
	 group=allele,
	 data=df,
	 auto.key=list(points = FALSE, rectangles = TRUE,space="top",size=2,cex=0.8),
	 stack=FALSE,
	 scales = scales,
	 ylab=gargs$ylab,
	 xlab=gargs$xlab,
	 box.ratio=2,
	 abbreviate=TRUE,
	 par.settings=parset,
	 main=amainVec
	)

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
#' 	r <- reads
#' 	seqlevels(r) <- "17"
#' 	covMatList <- coverageMatrixListFromGAL(BamList=r, strand="+")
#' 
#' @export coverageMatrixListFromGAL
coverageMatrixListFromGAL <- function(BamList,strand=NULL,ignore.empty.bam.row=TRUE){

	#If having common start and end points for all gviz track objects the matrix will start on the specific start regardless if there are reads in the bamList or not. 
	
	#TODO, to conveniently access data without loading into memory, the bam file should be read again by using the argument BamPath. 

	GAL <- BamList

	#Could be good with a check that the matrix is not longer than
	#CNTNAP2 - 2300000bp long which is the  longest gene
	#But will wait with that.
	pstrand=FALSE
	mstrand=FALSE
	if(!is.null(strand)){
		if(strand=="+"){pstrand=TRUE}
		else if(strand=="-"){mstrand=TRUE}
		else{stop("strand has to be '+' or '-' if not NULL\n")}
	}

	if(!length(seqlevels(GAL))==1){stop("can only be one seq level\n")}

	#get start and end before filtering on strand, will make things easier downstream.
	suppressWarnings(bamStart <- min(min(start(GAL))))
	suppressWarnings(bamEnd <- max(max(end(GAL))))
	bamWidth <- bamEnd-bamStart+1

	#if(is.null(start) | is.null(end)){
		start <- bamStart
		end <- bamEnd
		width <- bamWidth
	#}

	if(pstrand){GALp <- GAL[strand(GAL)=="+"]}
	if(mstrand){GALm <- GAL[strand(GAL)=="-"]}

	if(pstrand){matP <- matrix(0,ncol=(width),nrow=length(GAL))}
	if(mstrand){matM <- matrix(0,ncol=(width),nrow=length(GAL))}
	if(pstrand){rownames(matP) <- names(GAL)}
	if(mstrand){rownames(matM) <- names(GAL)}
	##################################################

	covVecFromGA <- function(GA){
			mcols(GA) <- NULL
			one <- unlist(grglist(GA))
			covRle <- coverage(one)[[1]]
			cov <- as.integer(window(covRle,start,end))
			cov
	}

	if(pstrand){for(i in 1:length(GALp)){matP[i,]<- covVecFromGA(GALp[[i]])}}
	if(mstrand){for(i in 1:length(GALm)){matM[i,]<- covVecFromGA(GALm[[i]])}}
	
	#make mat from matP or matM
	if(pstrand){mat <- matP }
	if(mstrand){mat <- matM }

	#store in a list
	if(!is.null(strand)){
		retList <- list(mat,start,end)
	}else{stop("strand must be present")}
	#set name on list
	if(!is.null(strand)){
		names(retList) <- c("mat","start","end")
	}else{stop("strand must be present")}
	
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
#' 	lst <- list(hungry="yes", thirsty="no")
#' 	implodeList(lst)
#'	#the check ls()
#'  ls()
#' @export implodeList
implodeList <- function(x){ 
	oname <- deparse(substitute(x))
	eval(parse(text=paste0("for(i in 1:length(",oname,")){assign(names(",oname,")[i],",oname,"[[i]])}")), 
		 parent.frame())
}



