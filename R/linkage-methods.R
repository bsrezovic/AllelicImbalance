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
#' @param rv riskVariant object with phase and 'ref'/'alt' allele information
#' @param region riskVariant object with phase and alternative allele information
#' @param settings riskVariant object with phase and alternative allele information
#' @param return.class 'vector' or 'matrix'
#' @param return.meta logical to return a list with metainformation
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
#' rv <- riskVariantFromGRangesAndPhaseArray(x=GRvariants, phase=p.ar)
#'
#' # in this example each and every snp in the ASEset defines a region
#' r1 <- granges(a)
#'
#' # in this example two overlapping subsets of snps in the ASEset defines the region
#' r2 <- split(granges(a)[c(1,2,2,3)],c(1,1,2,2))
#'
#' # link variant almlof (lva)
#' lva(a, rv, r1)
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
		if("threshold.distance" %in% names(settings)){
			distance <- settings[["threshold.distance"]]
		}else{
			distance <- 200000
		}

		#region summary
		rs <- regionSummary(x, region)
							

		#match riskVariant to rs granges
		hits <- findOverlaps(rv, rs$gr + distance)
		#stop if no overlap
		if(length(rv2)==0){stop("no overlap between rs and rv")}
		#if overlap, then subset hits
		rs2 <- rs$x[,,subjectHits(hits), drop=FALSE]
		rv2 <- rv[queryHits(hits), drop=FALSE]

		#make groups for regression based on (het hom het)
		grp <- .groupBasedOnPhaseAndAlleleCombination(phase(rv2,return.class="array")[,,c(1, 2), drop=FALSE])
		#call internal regression function	
		pvalues <- lva.internal(rs2, grp)

		#create an object with results
		#lva <- new("Lva", sset,
		#	meta = list()
		#)

		#create return object
		if(return.class=="vector"){
			pvalues
		}else if(return.class=="matrix"){
			if("ixn" %in% names(rs2)){
				rs2.names <- apply(rs$ixn[-nrow(rs$ixn),
								   subjectHits(hits)],2,paste,collapse="/")
				region.annotation <- rs2.names
			}else if("idx.names" %in% names(rs)){
				rs2.names <- rs$idx.names[subjectHits(hits)]
				region.annotation <- rs2.names
			}else{
				region.annotation <- "nosetname"
			}

			if(return.meta){
				list(mat=matrix(c(pvalues, rownames(rv2), region.annotation),ncol=3),
					 GRrv=rowRanges(rv2),
					 GRrs=rs$gr[subjectHits(hits)])
			}else{
				matrix(c(pvalues, rownames(rv2), region.annotation),ncol=3)
			}
		}
})


### -------------------------------------------------------------------------
### helpers for lva
###

#send in an array rows=SNPs, cols=sampels, 3dim=phase with maternal as el. 1 and paternal as el. 2)
#returns a matrix with three groups, het A/B =1, hom (A/A, B/B) =2 and het B/A =3
.groupBasedOnPhaseAndAlleleCombination <- function(ar){
		grp <- matrix(2, nrow=dim(ar)[2], ncol=dim(ar)[1])
		grp[(ar[,,1] == 1) & (ar[,,2] == 0)] <- 1                             	
		grp[(ar[,,1] == 0) & (ar[,,2] == 1)] <- 3
		grp
}

#' lva.internal
#' 
#' make an almlof regression for arrays (internal core function)
#' 
#' internal method that takes one array with results from regionSummary
#' and one matrix with group information for each risk SNP (based on phase)
#'
#' @name lva.internal
#' @rdname lva.internal
#' @aliases lva.internal,array-method
#' @docType methods
#' @param x regionSummary array phased for maternal allele
#' @param grp group 1-3 (1 for 0:0, 2 for 1:0 or 0:1, and 3 for 1:1)
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
#' lva.internal(rs, grp, 3)
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
		
		unlist(.lvaRegressionPvalue(x, grp, element))
})

### -------------------------------------------------------------------------
### helpers for lva.internal
###

# input ar array 1d=samples, 2d=variable, 3d=SNP
# input grp matrix with group 1 2 3. 2d=samples, 1d=SNP
# lapply over each snp and make regression over variable element based on grp
# output is a list with one result for each SNP
.lvaRegressionPvalue <- function(ar, grp, element){
	lapply(1:dim(ar)[3], function(i, y, x){
				summary(lm(y[, element, 1]~x[i, ]))$coefficients[2, 4]
		}, y=ar, x=grp)
}






