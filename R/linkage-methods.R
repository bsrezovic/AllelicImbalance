#'@include ASEset-class.R
NULL


#' lva
#' 
#' make an almlof regression for arrays
#' 
#' internal method that takes one array with results from regionSummary
#' and one matrix with group information for each risk SNP (based on phase)
#'
#' @name lva
#' @rdname lva
#' @aliases lva,array-method
#' @docType methods
#' @param x ASEset object with phase and 'ref'/'alt' allele information
#' @param rv RiskVariant object with phase and 'ref'/'alt' allele information
#' @param region RiskVariant object with phase and alternative allele information
#' @param settings RiskVariant object with phase and alternative allele information
#' @param return.class 'LinkVariantAlmlof' (more options in future)
#' @param verbose logical, if set TRUE, then function will be more talkative
#' @param ... arguments to forward to internal functions
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords phase
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
#' #init risk variants
#' p.ar <- phaseMatrix2Array(p)
#' rv <- RiskVariantFromGRangesAndPhaseArray(x=GRvariants, phase=p.ar)
#'
#' #colnames has to be samea and same order in ASEset and RiskVariant
#' colnames(a) <- colnames(rv)
#'
#' # in this example each and every snp in the ASEset defines a region
#' r1 <- granges(a)
#' 
#' #use GRangesList to merge and use regions defined by each element of the
#' #GRangesList
#' r1b <- GRangesList(r1)
#' r1c <- GRangesList(r1, r1)
#'
#' # in this example two overlapping subsets of snps in the ASEset defines the region
#' r2 <- split(granges(a)[c(1,2,2,3)],c(1,1,2,2))
#'
#' # link variant almlof (lva)
#' lva(a, rv, r1)
#' lva(a, rv, r1b)
#' lva(a, rv, r1c)
#' lva(a, rv, r2)
#' 
NULL

#' @rdname lva
#' @export
setGeneric("lva", function(x, ... 
	){
    standardGeneric("lva")
})

#' @rdname lva
#' @export
setMethod("lva", signature(x = "ASEset"),
		function(x, rv, region, settings=list(),
				 return.class="LinkVariantAlmlof",
				 verbose=FALSE, ...
	){

		#safety check
		if(any(!colnames(x) %in% colnames(rv)) | any(!colnames(rv) %in% colnames(x))) 
				stop("missmatch of colnames for x and rv")

		if("threshold.distance" %in% names(settings)){
			distance <- settings[["threshold.distance"]]
		}else{
			distance <- 200000
		}

		#region summary
		rs <- regionSummary(x, region)


		#match riskVariant to rs granges
		hits <- findOverlaps(rv, granges(rs) + distance)
		#stop if no overlap
		if(length(hits)==0){stop(paste("no rs and rv are close with distance: ",distance, sep=""))}
		#if overlap, then subset hits
		rs2 <- rs[subjectHits(hits)]
		rv2 <- rv[queryHits(hits),, drop=FALSE]
		#make groups for regression based on (het hom het)
		grp <- .groupBasedOnPhaseAndAlleleCombination(phase(rv2,return.class="array")[,,c(1, 2), drop=FALSE])
		plotGroups <- .lvaGroups(mcols(rv2)[["ref"]],mcols(rv2)[["alt"]], )
		#call internal regression function	
		mat <- lva.internal(assays(rs2)[["rs1"]], t(grp))

		#make txSNP specific lva test
		rs2 <- .addLva2ASEset(rs2, grp)

		#create return object
		if(return.class=="LinkVariantAlmlof"){
			sset <- SummarizedExperiment(
						assays = SimpleList(rs1=assays(rs2)[["rs1"]], lvagroup=grp), 
						colData = colData(rs2),
						rowRanges = granges(rs2))

			rownames(sset) <- rownames(rs2)
			mcols(sset)[["RiskVariantMeta"]] <- DataFrame(GR=granges(rv2), rsid=rownames(rv2))
			mcols(sset)[["RiskVariantMetaFull"]] <- rv2
			mcols(sset)[["LMCommonParam"]] <- DataFrame(mat, row.names=NULL)
			mcols(sset)[["LvaPlotGroups"]] <- DataFrame(plotGroups, row.names=NULL)

			#create an object with results
			new("LinkVariantAlmlof", sset,
				meta = list()
			)
		}
})


### -------------------------------------------------------------------------
### helpers for lva
###

#send in an array rows=SNPs, cols=sampels, 3dim=phase with maternal as el. 1 and paternal as el. 2)
#returns a matrix with three groups, het A/B =1, hom (A/A, B/B) =2 and het B/A =3
.groupBasedOnPhaseAndAlleleCombination <- function(ar){
		grp <- matrix(2, nrow=dim(ar)[1], ncol=dim(ar)[2])
		grp[(ar[,,1] == 1) & (ar[,,2] == 0)] <- 1                             	
		grp[(ar[,,1] == 0) & (ar[,,2] == 1)] <- 3
		grp
}


.lvaGroups <- function(ref, alt){
		fir <- paste(alt,"|", ref, sep="")
		sec <- paste(ref,"|", ref, " AND ",alt,"|", alt, sep="")
		thi <- paste(ref,"|", alt, sep="")

		matrix(c(fir, sec, thi), ncol=3)
}

.addLva2ASEset <- function(rs2, grp){
		lst <- mcols(rs2)[["ASEsetMeta"]][[1]]
		for(i in 1:length(lst)){
		  fr <- assays(lst[[i]])[["matfreq"]]
		  grp2 <- grp[i,]
		  lmcomparam <- .lvaRegressionReturnCommonParamMatrixTxSNPspecific(fr,grp2)
		  mcols(mcols(rs2)[["ASEsetMeta"]][[1]][[i]])[["lmcomparam"]] <- DataFrame(lmcomparam)
		}
		rs2
}


#' lva.internal
#' 
#' make an almlof regression for arrays (internal function)
#' 
#' internal method that takes one array with results from regionSummary
#' and one matrix with group information for each risk SNP (based on phase).
#' Input and output objects can change format slightly in future.
#'
#' @name lva.internal
#' @rdname lva.internal
#' @aliases lva.internal,array-method
#' @docType methods
#' @param x regionSummary array phased for maternal allele
#' @param grp group 1-3 (1 for 0:0, 2 for 1:0 or 0:1, and 3 for 1:1)
#' @param element which column in x contains the values to use with lm.
#' @param ... arguments to forward to internal functions
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords phase
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
#' # in this example two overlapping subsets of snps in the ASEset defines the region
#' region <- split(granges(a)[c(1,2,2,3)], c(1,1,2,2))
#' rs <- regionSummary(a, region, return.class="array", return.meta=FALSE)
#'
#' # use  (change to generated riskSNP phase later)
#' phs <- array(c(phase(a,return.class="array")[1,,c(1, 2)], 
#'				 phase(a,return.class="array")[2,,c(1, 2)]), dim=c(20,2,2))
#' grp <- matrix(2, nrow=dim(phs)[1], ncol=dim(phs)[2])		 
#' grp[(phs[,,1] == 0) & (phs[,,2] == 0)] <- 1
#' grp[(phs[,,1] == 1) & (phs[,,2] == 1)] <- 3

#' #only use mean.fr at the moment, which is col 3
#' lva.internal(assays(rs)[["rs1"]], grp, 3)
#' 
NULL

#' @rdname lva.internal
#' @export
setGeneric("lva.internal", function(x, ... 
	){
    standardGeneric("lva.internal")
})

#' @rdname lva.internal
#' @export
setMethod("lva.internal", signature(x = "array"),
		function(x, grp, element=3, ...
	){
		
		#unlist(.lvaRegressionPvalue(x, grp, element))
		.lvaRegressionReturnCommonParamMatrix(x, grp, element)

})

### -------------------------------------------------------------------------
### helpers for lva.internal
###

# input ar array 1d=SNP, 2d=samples, 3d=variable
# input grp matrix with group 1 2 3. 2d=samples, 1d=SNP
# lapply over each snp and make regression over variable element based on grp
# output is a list with one result for each SNP
.lvaRegressionPvalue <- function(ar, grp, element){
	lapply(1:dim(ar)[1], function(i, y, x){
				summary(lm(y[i, ,element]~x[, i]))$coefficients[2, 4]
		}, y=ar, x=grp)
}



.lvaRegressionReturnCommonParamMatrix <- function(ar, grp, element){

	mat <- matrix(NA, ncol=dim(ar)[3], nrow=nrow(ar))
	nocalc <- apply(ar[,,3, drop=FALSE], 1, function(x){sum(!(is.na(x)))==0})

	#only make regression if there is at least one row possible to compute
	if(any(!nocalc)){
		mat[!nocalc,] <- t(sapply(which(!nocalc), function(i, y, x){
						mat2 <- matrix(NA, ncol=2, nrow=4)
						s <-summary(lm(y[i, ,element]~x[, i]))$coefficients
						mat2[,1:nrow(s)] <- s
						c(mat2)
					}, y=ar[!nocalc,,,drop=FALSE], x=grp[,!nocalc,drop=FALSE]))
	}
	colnames(mat) <- c("est1","est2","stderr1","stderr2","tvalue1","tvalue2","pvalue1","pvalue2")
	mat
}

.lvaRegressionReturnCommonParamMatrixTxSNPspecific <- function(fr, grp){
	fr2 <- t(fr)
	grp2 <- grp
	mat <- matrix(NA, ncol=8, nrow=ncol(fr2))
	nocalc <- apply(fr2[,, drop=FALSE], 2, function(x){sum(!(is.na(x)))==0})
		#y <- fr2[!nocalc,,drop=FALSE]
		#x <- grp[,!nocalc,drop=FALSE]
	    x <- grp2
		for(i in 1:ncol(fr2)){
			y <- fr2[,i]
			mat2 <- matrix(NA, ncol=2, nrow=4)
			if(!(nocalc[i])){
			  s <-summary(lm(y~x))$coefficients
			  mat2[,1:nrow(s)] <- s
			  mat[i,] <- c(mat2)
			}
		}
	colnames(mat) <- c("est1","est2","stderr1","stderr2","tvalue1","tvalue2","pvalue1","pvalue2")

	#only make regression if there is at least one row possible to compute
	#if(any(!nocalc)){
	#	mat[!nocalc,] <- t(sapply(which(!nocalc), function(i, y, x){
	#					mat2 <- matrix(NA, ncol=2, nrow=4)
	#					s <-summary(lm(y[i, ]~x[, i]))$coefficients
	#					mat2[,1:nrow(s)] <- s
	#					c(mat2)
	#				}, y=fr2[!nocalc,,drop=FALSE], x=grp[,!nocalc,drop=FALSE]))
	#}

	mat
}



