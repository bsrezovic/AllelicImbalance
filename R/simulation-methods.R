
##the draft of the function giving us a fastq file for all samples to detect mapbias
#
##' simulate reads from ASEset
##' 
##' simulate reads
##' 
##' to compensate for reference genome mapping bias
##' 
##' @param x ASEset
##' @param reads number of reads to create
##' @param refPath reference genome path
##' @param read.length length of reads
##' @param PE paired-end reads, TRUE or FALSE
##' @param insert.size Nr bases separating the two reads, only if PE=TRUE 
##' @param insert.sd standard deviation for the insert.size
##' @author Jesper R. Gadin
##' @keywords infer
##' @examples
##' 
##' data(ASEset)
##' g <- makeMapBiasReads(ASEset)
##' 
##' @export makeMapBiasReads
#
#makeMapBiasReads <- function(x, reads=100, refPath, read.length=100, PE=FALSE, insert.size=-30, insert.sd=1){
#
#	#tmp ref file
#
#
#	#the output will be a GRanges or GrangesList  with reads
#
#	# x is of ASEset class
#	# reads are the numeber of simulated reads for each variant
#	# ref is the reference genome used
#
#	#require
#	require(Rsamtools)
#
#	#no junctions(if junctions we need txdb)
#	searchArea.ref <- rowData(x) + read.length-1
#
#	#make FaFile
#	fl <- FaFile("hg19.fa")
#
#	#check if the index file is present otherwise tell the user to use the indexFa(FaFile("pathToReference")) command
#	if(list.files("placeholder", pattern=paste("*",ref,".fai",sep=""))){
#		indexFa(FaFile("pathToReference"))
#	}
#	#IMPORTANT! The  index command only needs to be executed once
#	#indexFa(fl) #creates a new file as index in the same directory but with extension *.fai 
#
#	#open,scan,close file
#	open(fl)
#	ref <- scanFa(fl,param=rowData(x))
#	seq <- scanFa(fl,param=searchArea.ref)
#	close(fl)
#	
#
#
#	#do it for one sample at the time
#	if(PE){
#		#for each snp find all neighbouring snps to that one (all within read size)
#		#and which of these are different from the ref genome?
#		for(i in 1:nrow(x)){
#			hits <- findOverlaps(rowData(x)[i]+read.length, rowData(x))				
#			#here we need the functionality to store the genotype information for each 
#		       	#snp. So while doing that, I might already plan to include all genotypes and their ranges in a slot that. It would be some kind of transposed colData table. In the end the genotypes for RNA-seqenced allele we can store in an assay	
#			(x)
#		}
#	}
#}



##construct the paried data table
##known phases
#phaseA <- c("TGACT")
#phaseB <- c("CTTGA")
#
##known pairs, and known phase
#pairA <- c("TG","GA","AC","CT")
#pairB <- c("CT","TT","TG","GA")
#
##known pairs, but unknown phase
##X and Y can be seen as temporary phases)
##that the pairs have randomly been located to
##but kepping the location information for each pair
#pairX <- c("TG","TT","AC","CT")
#pairY <- c("CT","GA","TG","GA")
#
##split pair
#pairX.split <- unlist(strsplit(pairX,split="",fixed=TRUE))
#pairY.split <- unlist(strsplit(pairY,split="",fixed=TRUE))
#
##first position in pair
#px1 <- pairX.split[seq.int(1L,length(pairX.split),2L)]
#py1 <- pairY.split[seq.int(1L,length(pairY.split),2L)]
##remove first pair
#px1 <- px1[-1]
#py1 <- py1[-1]
#
##second position in pair
#px2 <- pairX.split[seq.int(2L,length(pairX.split),2L)]
#py2 <- pairY.split[seq.int(2L,length(pairY.split),2L)]
##remove first pair
#px2 <- px2[-c(length(px2))]
#py2 <- py2[-c(length(py2))]
#
##split pairs tables in array (third dimension is for each pair)
#snpPairs <- length(px1)
#ar <- array(c(px2,py2,px2,py2,px1,py1,py1,px1),dim=c(snpPairs,4,2))
#ar <- aperm(ar,c(2,3,1))
#
##check which pairs belong to eachother
#w <- which(ar[,1,] ==ar[,2,]
#)
#
##take out odd numbers which reflects the X strand and then the Y strand (even numbers)
#X <- w[w %% 2!=0]
#Y <- w[w %% 2==0]
#
##move all numbers down to baseline 1-4
#X2 <- X-seq(0,length(X)*3,4)
#Y2 <- X-seq(0,length(Y)*3,4)
#
##compensate for the disorder when one pair was not one the X template
#Rl <- Rle(X2)
#v <- X2
#v[v == 3 ] <- unlist(lapply(runLength(Rl)[runValue(Rl) == 3 ], 
#        function(x) {
#            c(3, rep(1, x - 1))
#        }))
#X2 <- v
#
##phase genotype 1 is same as X template 3 is same as Y Template (X always starts with 1 or TRUE)
#X.TF <- c(TRUE,X2==1)
#X.Ph <- pairX
#X.Ph[!X.TF] <- pairY[!X.TF]
#
##validate to the known phase
#identical(X.Ph,pairA)
#
