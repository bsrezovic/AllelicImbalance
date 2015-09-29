#'@include ASEset-class.R
NULL

#' regionSummary
#' 
#' Gives a summary of AI-consistency for a transcript
#'
#' From a given set of e.g. transcripts exon ranges the function will return
#' a summary for the sum of all exons. Phase information, reference and alternative
#' allele is required.
#'
#' A limitation comes to the strand-specificness. At the moment it is not possible
#' to call over more than one strand type using the strands in region. This will be
#' improved before going to release. 
#'
#' to calculate the direction and binomial p-values of AI the mapbias stored in the
#' ASEset is used. see '?mapBias'.
#'
#'
#' @name regionSummary
#' @rdname regionSummary
#' @aliases regionSummary,numeric-method
#' @docType methods
#' @param x ASEset object
#' @param region to summmarize over, the object can be a GRanges, GRangesList 
#' @param strand can be "+", "-" or "*"
#' @param return.class "array" or "list".
#' @param ... arguments to forward to internal functions
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords summary
#' @examples
#' 
#' data(ASEset) 
#' a <- ASEset
#' # Add phase
#' set.seed(1)
#' p1 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
#' p2 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
#' p <- matrix(paste(p1,sample(c("|","|","/"), size=nrow(a)*ncol(a), replace=TRUE), p2, sep=""),
#' 	nrow=nrow(a), ncol(a))
#' 
#' phase(a) <- p
#' 
#' #add alternative allele information
#' mcols(a)[["alt"]] <- inferAltAllele(a)
#'
#' # in this example each and all snps in the ASEset defines the region
#' region <- granges(a)
#' t <- regionSummary(a, region)
#'
#' # in this example two overlapping subsets of snps in the ASEset defines the region
#' region <- split(granges(a)[c(1,2,2,3)],c(1,1,2,2))
#' t <- regionSummary(a, region)
#'
NULL

#' @rdname regionSummary
#' @export
setGeneric("regionSummary", function(x, ... ){
    standardGeneric("regionSummary")
})


#' @rdname regionSummary
#' @export
setMethod("regionSummary", signature("ASEset"),
		function(x, region, strand="*", 
				 return.class="RegionSummary", ...
	){

		#checks are needed for
		#1 phase
		#2 alt ALT
		#3 ref REF

		if(class(region)=="GRanges" | class(region)=="GRangesList"){
			reg <- .unlistGRangesListAndIndex(region)
		}else if(class(region)=="list"){cat("list version not implemented")}

		#make overlap and subset for ASEset
		x <- .selectRegionAndTransferIndexToASEset(x, reg)
		#make regional granges for she summaries
		gr <- .makeRegionGRangesFromASEsetWithRegionIndex(x)
		#fraction with maternal phase used as numerator
		fr <- fraction(x, strand=strand, top.fraction.criteria="phase")
		#which snp and sample pairs are heterozygotes
		fr.het.filt <- hetFilt(x)
		#filter so all homozygotes become NAs 
		fr[!t(fr.het.filt)] <- NA
		#maternal phase map bias
		mb <- t(.maternalPhaseMapBias(mb=mapBias(x, return.class="array"), ma=maternalAllele(x), va=x@variants))
		#calculate delta
		fr.d <- .deltaFromFractionMatrixAndMapBias(fr, mb)
		#pick our important variables
		idx <-  mcols(x)[["regionIndex"]][[1]]
		#Summarize the statistics
		hets <- .regionStatisticsFromMatrixAndIndex(t(fr.het.filt), idx, sum)
		homs <- .regionStatisticsFromMatrixAndIndex(t(!fr.het.filt), idx, sum)
		mean.delta <- .regionStatisticsFromMatrixAndIndex(fr.d, idx, mean)
		sd.delta <- .regionStatisticsFromMatrixAndIndex(fr.d, idx, mean)
		mean.fr <- .regionStatisticsFromMatrixAndIndex(fr, idx, mean)
		sd.fr <- .regionStatisticsFromMatrixAndIndex(fr, idx, mean)
		ai.dir <- .aiRegionDirectionFromMaternalPhaseMapBiasAndIndex(fr, mb, idx)

		#For Meta infor, split ASEset and put into a DataFrame with same nrows as return object 
		xdf <- .makeMetaASEsetDataFrameForRegionSummary(x, idx)

		#make array
		ar <- aperm(array(c(hets, homs, mean.fr, sd.fr, mean.delta, sd.delta, ai.dir),
			  dim=c(ncol(x),length(unique(idx)), 8),), c(2,1,3))
		
		dimnames(ar) <- list(1:length(unique(idx)), colnames(x), 
			  c("hets","homs","mean.fr","sd.fr","mean.delta","sd.delta","ai.up","ai.down"))

		#create RegionSummary object
		sset <- SummarizedExperiment(
					assays = SimpleList(rs1=ar), 
					colData = colData(x),
					rowRanges = gr)

		rownames(sset) <- names(gr)

		#add ASEset meta objects
		mcols(sset)[["ASEsetMeta"]] <- xdf

		#valid
		#validObject(.Object)

		#Return object
		new("RegionSummary", sset,
			meta = list(),
			sumnames = dimnames(ar)[[3]]
		)
})


### -------------------------------------------------------------------------
### helpers for regionSummary
###

# fr is the value for when maternal phase is the numerator in the freq. ekv.
# output: dim3 first element is up and second element is down
.aiRegionDirectionFromMaternalPhaseMapBiasAndIndex<- function(fr, mb, idx){

		fr.up.filt <- fr > mb
		fr.down.filt <- fr < mb

		mode(fr.up.filt) <- "integer"
		mode(fr.down.filt) <- "integer"

		array(c(
			 t(rowsum(t(fr.up.filt),idx, na.rm=TRUE)), 
			 t(rowsum(t(fr.down.filt),idx, na.rm=TRUE))
		 ), dim=c(nrow(fr), length(unique(idx)), 2))
}

.maternalPhaseMapBias<- function(mb, ma, va){
		#for loop (for now)
		vmat <- matrix(va, byrow=TRUE, ncol=length(va), nrow=ncol(mb))
		mbias.values <- matrix(NA, ncol=ncol(mb), nrow=nrow(mb))
		for (i in 1:nrow(mb)){
			it.mbias <- mb[i,,]
			it.mallele <- ma[i,]
			mat.mallele <- matrix(it.mallele, ncol=length(va), nrow=nrow(vmat))
			tf <- mat.mallele == vmat
			mbias.values[i,] <- it.mbias[tf]
		}
		mbias.values
}

#input: ar has dim(samples, SNPs)
.regionStatisticsFromMatrixAndIndex<- function(mat, idx, fun){
			t(tapply(t(mat), list(idx[row(t(mat))], col(t(mat))), fun, na.rm=TRUE)) 
		
}

.deltaFromFractionMatrixAndMapBias<- function(ar, mb){
			abs(ar - mb)		
}

# input, ASEset and unlisted region, (region is a GRanges object, with index information)
# output an indexed and sorted(lowest idx to highest) ASEset overlapping the region
.selectRegionAndTransferIndexToASEset <- function(a, reg){

		hits <- findOverlaps(a, reg)
		a <- a[queryHits(hits)]
		mcols(a)[["ASEsetIndex"]] <- queryHits(hits)
		mcols(a)[["regionIndex"]] <- mcols(reg[subjectHits(hits)])[["regionIndex"]]
		mcols(a)[["regionIndexName"]] <- mcols(reg[subjectHits(hits)])[["regionIndexName"]]
		a <- a[sort(mcols(a)[["regionIndex"]][[1]], index.return=TRUE)$ix]
		a
}

#sorted and indexed ASEset
.unlistGRangesListAndIndex <- function(grl){
		idx <- togroup(PartitioningByWidth(elementLengths(grl)))
		if(!is.null(names(grl))){
			idn <- names(grl)[idx]
		}else{
			idn <- "nameless"
		}
		new <- unlist(grl)
		mcols(new)[["regionIndex"]] <- DataFrame(lvl1=idx)
		mcols(new)[["regionIndexName"]] <- DataFrame(lvl1=idn)
		new
}

#return granges for each bin (index has to be sorted)
.makeRegionGRangesFromASEsetWithRegionIndex <- function(x){
	idx <- mcols(x)[["regionIndex"]][[1]]
	idn <- mcols(x)[["regionIndexName"]][[1]]
	spl <- split(granges(x), idx)
	il <- IntegerList(lapply(spl,function(x){mcols(x)[["ASEsetIndex"]]}))
	names(il) <- NULL
	gr <- unlist(reduce(spl, min.gapwidth=100000000))
	mcols(gr) <- NULL
	mcols(gr)[["ASEsetIndex"]] <- il
	mcols(gr)[["regionIndex"]] <- DataFrame(lvl1=as.integer(unique(idx)))
	mcols(gr)[["regionIndexName"]] <- DataFrame(lvl1=as.character(unique(idn)))
	gr
}

#make list of the ASEset components relative to region
.makeMetaASEsetDataFrameForRegionSummary <- function(as, idx){
		DataFrame(ASEsetList=split(as, idx))
}



