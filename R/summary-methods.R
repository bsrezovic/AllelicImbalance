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
#' @param return.meta logical if to return a list with additional metadata
#' @param drop logical that for TRUE drops the third dimentsion if there is only
#' one element.
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
#' # in this example all snps in the ASEset defines the region
#' region <- granges(a)
#' t <- regionSummary(a, region)
#'
#' # in this example two overlapping subsets of snps in the ASEset defines the region
#' region <- split(granges(a)[c(1,2,2,3)],c(1,1,2,2))
#' t <- regionSummary(a, region)
#'
#' # use a multilevel list as input (output will keep the list dimensions)
#' region <- split(granges(a)[c(1,2,2,3)],c(1,1,2,2))
#' names(region) <- c("introns", "exons")
#' region <- list(g1=list(tx1=region, tx2=region), g2=list(tx1=region, tx2=region, tx3=region))
#' t <- regionSummary(a, region)
NULL

#' @rdname regionSummary
#' @export
setGeneric("regionSummary", function(x, ... ){
    standardGeneric("regionSummary")
})


#' @rdname regionSummary
#' @export
setMethod("regionSummary", signature("ASEset"),
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

