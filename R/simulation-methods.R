
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
