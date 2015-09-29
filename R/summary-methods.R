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
<<<<<<< HEAD
#' improved before going to release. 
#'
#' to calculate the direction and binomial p-values of AI the mapbias stored in the
#' ASEset is used. see '?mapBias'.
#'
=======
#' improved before going to release.
>>>>>>> chopped out linkage and summary methods R files and chopped out a first batch of units and tested them
#'
#' @name regionSummary
#' @rdname regionSummary
#' @aliases regionSummary,numeric-method
#' @docType methods
#' @param x ASEset object
#' @param region to summmarize over, the object can be a GRanges, GRangesList or 
#' list containing a GRangesList in the root
#' @param strand can be "+", "-" or "*"
#' @param threshold.pvalue used in filter when to count AI as significant
#' @param return.class "array" or "list".
<<<<<<< HEAD
=======
#' @param return.meta logical if to return a list with additional metadata
#' @param drop logical that for TRUE drops the third dimentsion if there is only
#' one element.
>>>>>>> chopped out linkage and summary methods R files and chopped out a first batch of units and tested them
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
<<<<<<< HEAD
#' # in this example each and all snps in the ASEset defines the region
=======
#' # in this example all snps in the ASEset defines the region
>>>>>>> chopped out linkage and summary methods R files and chopped out a first batch of units and tested them
#' region <- granges(a)
#' t <- regionSummary(a, region)
#'
#' # in this example two overlapping subsets of snps in the ASEset defines the region
#' region <- split(granges(a)[c(1,2,2,3)],c(1,1,2,2))
#' t <- regionSummary(a, region)
#'
<<<<<<< HEAD
#' 
=======
#' # use a multilevel list as input (output will keep the list dimensions)
#' region <- split(granges(a)[c(1,2,2,3)],c(1,1,2,2))
#' names(region) <- c("introns", "exons")
#' region <- list(g1=list(tx1=region, tx2=region), g2=list(tx1=region, tx2=region, tx3=region))
#' t <- regionSummary(a, region)
>>>>>>> chopped out linkage and summary methods R files and chopped out a first batch of units and tested them
NULL

#' @rdname regionSummary
#' @export
setGeneric("regionSummary", function(x, ... ){
    standardGeneric("regionSummary")
})


#' @rdname regionSummary
#' @export
setMethod("regionSummary", signature("ASEset"),
<<<<<<< HEAD
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
		#pick our important variables
		idx <-  mcols(x)[["regionIndex"]][[1]]
		#fraction with maternal phase used as numerator
		fr <- fraction(x, strand=strand, top.fraction.criteria="phase")
		#which snp and sample pairs are heterozygotes
		fr.het.filt <- hetFilt(x)
		#filter so all homozygotes become NAs 
		fr[!t(fr.het.filt)] <- NA
		#maternal phase map bias
		print(dim(x))
		mb <- t(.maternalPhaseMapBias(mb=mapBias(x, return.class="array"), ma=maternalAllele(x), va=x@variants))
		#calculate delta
		print(dim(fr))
		print(dim(mb))
		fr.d <- .deltaFromFractionMatrixAndMapBias(fr, mb)
		
		#Summarize the statistics
		hets <- .regionStatisticsFromMatrixAndIndex(t(fr.het.filt), idx, sum)
		homs <- .regionStatisticsFromMatrixAndIndex(t(!fr.het.filt), idx, sum)
		mean.delta <- .regionStatisticsFromMatrixAndIndex(fr.d, idx, mean)
		sd.delta <- .regionStatisticsFromMatrixAndIndex(fr.d, idx, mean)
		mean.fr <- .regionStatisticsFromMatrixAndIndex(fr, idx, mean)
		sd.fr <- .regionStatisticsFromMatrixAndIndex(fr, idx, mean)
		ai.dir <- .aiRegionDirectionFromMaternalPhaseMapBiasAndIndex(fr, mb, idx)

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
		idn <- names(grl)[idx]
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

=======
		function(x, region, strand="*", threshold.pvalue=0.05, 
				 return.class="list", return.meta=TRUE, drop=TRUE, ...
	){

		#needs alternative allele
		if(!"alt" %in% colnames(mcols(x))){
			stop("function needs mcols(x)[['alt']] to be set")				
		}

		#need information of which are heterozygotes and homozygotes
		if(is.null(genotype(x))){
			genotype(x) <- inferGenotypes(x, return.allele.allowed="bi")
		}

		populate.list <- FALSE

		#check class of region
		if(class(region)=="GRanges"){
			idx <- rep(1,length(region))
			ar.dim3 <- 1
			ar.dim3.names <- "nameless"
		}else if(class(region)=="GRangesList")	{
			idx <- togroup(PartitioningByWidth(elementLengths(region)))
			ar.dim3 <- length(region)
			ar.dim3.names <- names(region)
		}else if(class(region)=="list") {
			idx.mat <- .multiUnlist.index(region)	
			idx.mat.names <- .multiUnlist.index.names(region)	
			rownames(idx.mat) <- paste("lvl", (nrow(idx.mat)+1):2, sep="")
			rownames(idx.mat.names) <- paste("lvl", (nrow(idx.mat.names)+1):2, sep="")
			region <- .multiUnlist(region)		
			populate.list <- TRUE

			idx <- togroup(PartitioningByWidth(elementLengths(region)))
			ar.dim3 <- length(region)
			ar.dim3.names <- 1:ar.dim3
		}

		#make overlap and subset based on gr
		hits <- findOverlaps(x, region)
		x <- x[queryHits(hits)[sort(subjectHits(hits), index.return=TRUE)$ix], ]

		if(return.meta){
			#return granges for each bin
			gr <- unlist(reduce(region, min.gapwidth=100000000))
		}

		fr.f <- fraction(x, strand=strand, top.fraction.criteria="phase")
		fr.het.filt <- hetFilt(x)

		#filter on p-value if threshold <1
		if(threshold.pvalue >= 1){
			fr.f[!t(fr.het.filt)] <- NA
		}else if(threshold.pvalue < 1 & threshold.pvalue > 0){
			pv <- binom.test(x,strand)
			fr.f[!(pv < threshold.pvalue) | !t(fr.het.filt)] <- NA
		}else{
			stop("threshold.pvalue must be a value > 0")
		}

		#maternal allele
		mallele <- maternalAllele(x)
		mbias <- mapBias(x, return.class="array")

		#for loop (for now)
		vmat <- matrix(x@variants, byrow=TRUE,ncol=length(x@variants),nrow=ncol(x))
		mbias.values <- matrix(NA, ncol=ncol(x), nrow=nrow(x))
		for (i in 1:nrow(mbias)){
			it.mbias <- mbias[i,,]
			it.mallele <- mallele[i,]
			mat.mallele <- matrix(it.mallele, ncol=length(x@variants), nrow=nrow(vmat))
			tf <- mat.mallele == vmat
			mbias.values[i,] <- it.mbias[tf]
		}

		fr.up.filt <- fr.f > t(mbias.values)
		fr.down.filt <- fr.f < t(mbias.values)

		mode(fr.up.filt) <- "integer"
		mode(fr.down.filt) <- "integer"
		
		#calc data
		ai.down <- t(rowsum(t(fr.down.filt),idx, na.rm=TRUE))
		ai.up <- t(rowsum(t(fr.up.filt),idx, na.rm=TRUE))
		fr.d <- abs(fr.f - t(mbias.values))
		hets <- t(tapply(fr.het.filt, list(idx[row(fr.het.filt)], col(fr.het.filt)), sum, na.rm=TRUE))
		homs <- t(tapply(!fr.het.filt, list(idx[row(!fr.het.filt)], col(!fr.het.filt)), sum, na.rm=TRUE))
		mean.fr <- t(tapply(t(fr.f), list(idx[row(t(fr.f))], col(t(fr.f))), mean, na.rm=TRUE))
		sd.fr <- t(tapply(t(fr.f), list(idx[row(t(fr.f))], col(t(fr.f))), sd, na.rm=TRUE))
		mean.delta <- t(tapply(t(fr.d), list(idx[row(t(fr.d))], col(t(fr.d))), mean, na.rm=TRUE))
		sd.delta <- t(tapply(t(fr.d), list(idx[row(t(fr.d))], col(t(fr.d))), sd, na.rm=TRUE))

		#make array
		ar <- aperm(array(c(hets, homs, mean.fr, sd.fr, mean.delta, sd.delta, ai.up, ai.down),
			  dim=c(ncol(x),ar.dim3,8),
			  dimnames=list(
						colnames(x),
						ar.dim3.names,
						c("hets",
						  "homs",
						  "mean.fr",
						  "sd.fr",
						  "mean.delta",
						  "sd.delta",
						  "ai.up",
						  "ai.down")
						 )
		), c(1,3,2))

		if(return.class=="array"){
			if(dim(ar)[3]==1 & drop){
				ar[,,1]
			}else{
				#check if index should be returned (recommended when wrapping GRangesList in lists )
				if(return.meta){
					if(populate.list){
						list(x=ar,ix=idx.mat,ixn=idx.mat.names, gr=gr)
					}else{
						list(x=ar, gr=gr, idx.names=ar.dim3.names)
					}
				}else{
					ar
				}
			}
		}else if(return.class=="list"){
			if(populate.list){
				lst <- .region.list.populate(ar, idx.mat[-nrow(idx.mat),], idx.mat.names[-nrow(idx.mat.names),])
				lst
			}else{

				lst <- lapply(seq(dim(ar)[3]), function(x) ar[ , , x])
				names(lst) <- ar.dim3.names

				if(length(lst)==1 & drop){
					lst[[1]]
				}else{
					lst
				}
			}
		}
})
>>>>>>> chopped out linkage and summary methods R files and chopped out a first batch of units and tested them

