#' phaseMatrix2Array
#' 
#' used to convert the phase from the visually friendly matrix to array.
#' 
#' A more effectice way of store the phase data in the ASEset object
#'
#' @name phaseMatrix2Array
#' @rdname phaseMatrix2Array
#' @aliases phaseMatrix2Array,matrix-method
#' @docType methods
#' @param x matrix see examples 
#' @param ... arguments to forward to internal functions
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords phase
#' @examples
#' 
#' #load data
#' data(ASEset) 
#' a <- ASEset
#'
#' #example phase matrix 
#' p1 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
#' p2 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
#' p <- matrix(paste(p1,sample(c("|","|","/"), size=nrow(a)*ncol(a), replace=TRUE), p2, sep=""),
#' 	nrow=nrow(a), ncol(a))
#'
#' ar <- phaseMatrix2Array(p)
#' 
NULL

#' @rdname phaseMatrix2Array
#' @export
setGeneric("phaseMatrix2Array", function(x, ... 
	){
    standardGeneric("phaseMatrix2Array")
})

#' @rdname phaseMatrix2Array
#' @export
setMethod("phaseMatrix2Array", signature(x = "matrix"),
		function(x, ...
	){

		psplit <- strsplit(x, split="")
		upsplit <- unlist(psplit)
		mat <- as.integer(upsplit[seq(1, length(upsplit), by=3)])
		pat <- as.integer(upsplit[seq(3, length(upsplit), by=3)])
		phased <- as.integer(upsplit[seq(2, length(upsplit), by=3)]=="|")
	
		array(c(mat,pat,phased), dim=c(nrow(x), ncol(x), 3))

})

#' phaseArray2Matrix
#' 
#' used to convert the phase from the visually friendly matrix to array.
#' 
#' A more effectice way of store the phase data in the ASEset object
#'
#' @name phaseArray2Matrix
#' @rdname phaseArray2Matrix
#' @aliases phaseArray2Matrix,array-method
#' @docType methods
#' @param x array see examples 
#' @param ... arguments to forward to internal functions
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords phase
#' @examples
#' 
#' #load data
#' data(ASEset) 
#' a <- ASEset
#'
#' #example phase matrix 
#' p1 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
#' p2 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
#' p <- matrix(paste(p1,sample(c("|","|","/"), size=nrow(a)*ncol(a), replace=TRUE), p2, sep=""),
#' 	nrow=nrow(a), ncol(a))
#' 
#' ar <- phaseMatrix2Array(p)
#'
#' #Convert back 
#' mat <- phaseArray2Matrix(ar)
#'
NULL

#' @rdname phaseArray2Matrix
#' @export
setGeneric("phaseArray2Matrix", function(x, ... 
	){
    standardGeneric("phaseArray2Matrix")
})

#' @rdname phaseArray2Matrix
#' @export
setMethod("phaseArray2Matrix", signature(x = "array"),
		function(x, ...
	){

		phased <- x[,,3]
		phased[phased==1] <- "|"
		phased[phased==0] <- "/"

		matrix(paste(x[,,1], phased, x[,,2], sep=""),
			nrow=nrow(x), ncol(x))
	
})

#' Plot Dataframe
#' 
#' Summarizes information to ease creating plots
#' 
#' Main purpose is to reduce the amount of overall code and ease maintenance. 
#'
#' top.allele.criteria can take three options, maxcount, ref and phase. The top
#' allele will be every second row in the data frame, with start from row 2. 
#' The maxcount argument will put the allele with most reads on top of the 
#' bivariate fraction. Similarly the ref argument will put always the reference
#' allele on top. The phase arguments puts the maternal phase always on top.
#' The top.allele.criteria for the ref or phase arguments requires that both ref
#' and alt is set in mcols(ASEset).
#' 
#'
#' @name fractionPlotDf
#' @rdname fractionPlotDf
#' @aliases fractionPlotDf,ASEset-method
#' @docType methods
#' @param x ASEset
#' @param snp rownames identifier for ASEset or row number
#' @param strand '+', '-' or '*'
#' @param top.allele.criteria 'maxcount', 'ref' or 'phase'
#' @param ... arguments to forward to internal functions
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords phase plotDf
#' @examples
#' 
#' #test on example ASEset 
#' data(ASEset)
#' a <- ASEset
#' df <- fractionPlotDf(a, 1, strand="+")
#'
#' @exportMethod fractionPlotDf
NULL

#' @rdname fractionPlotDf
#' @export
setGeneric("fractionPlotDf", function(x, snp,  strand="*", top.allele.criteria="maxcount", ... 
	){
    standardGeneric("fractionPlotDf")
})

#' @rdname fractionPlotDf
#' @export
setMethod("fractionPlotDf", signature(x = "ASEset"),
		function(x, snp,  strand="*", top.allele.criteria="maxcount", ...
	){

	##top.allele.criteria
	# maxcount
	# ref
	# phase
	
	if(class(snp)=="integer"){
		snprow <- snp
	}else if(class(snp)=="character"){
		snprow <- which(rownames(x) %in% snp  )
	}else{
		# If the class of snp is not approprite there will be an informative error message from
		# the base subset method.
		snprow <- snp
	}

	afraction <- fraction(x[snprow], strand = strand, top.allele.criteria=top.allele.criteria)

	values <- as.vector(t(matrix(as.numeric(afraction), ncol=2, nrow=ncol(x))))
	#values[seq(1,ncol(x)*2,by=2)+1] <- 1 - values[seq(1,ncol(x)*2,by=2)+1]

	#The first value will be the bottom allele (second allele the top allele)
	values[seq(1,ncol(x)*2,by=2)] <- 1 - values[seq(1,ncol(x)*2,by=2)]
	samples <- as.vector(t(matrix(rownames(afraction), ncol=2, nrow=ncol(x))))

	if(top.allele.criteria=="phase"){

		if(!sum(c("ref", "alt") %in% names(mcols(x))) == 2 ){
			stop("ref and alt has to be set in mcols to use phase option")
		}

		mat <- phase(x[snprow],return.class="array")[,,1]
		mat2 <- matrix(c(rep(mcols(x[snprow])[,"ref"],ncol(x)),rep(mcols(x[snprow])[,"alt"], ncol(x))), ncol=2, nrow=ncol(x))
		mat2[mat==0,1] <- mcols(x[snprow])[,"alt"]
		mat2[mat==0,2] <- mcols(x[snprow])[,"ref"]

		alleles <- as.vector(t(mat2))
		phase <- rep(c("paternal","maternal"), ncol(x))

	}else if(top.allele.criteria=="ref"){
		
		if(!sum(c("ref", "alt") %in% names(mcols(x))) == 2 ){
			stop("ref and alt has to be set in mcols to use phase option")
		}

		alleles <- as.vector(t(matrix(c(rep(mcols(x[snprow])[,"alt"],ncol(x)),rep(mcols(x[snprow])[,"ref"], ncol(x))), ncol=2, nrow=ncol(x))))
		
	}else if(top.allele.criteria=="maxcount"){

		#arank(x[snprow],return.class="matrix")[c(1,2)]
		alleles <- as.vector(t(matrix(arank(x[snprow], strand=strand, return.class="matrix")[c(2,1)], ncol=2, nrow=ncol(x), byrow=TRUE)))
		#alleles <- as.vector(t(matrix(arank(x[snprow], strand=strand, return.class="matrix")[c(1,2)], ncol=2, nrow=ncol(x), byrow=TRUE)))

	}

    TFna <- is.na(values)
    values[TFna] <- 0  # 0.5 + 0.5 -> 1
    na <- rep("no", length(values))
    na[TFna] <- "yes"

    df <- data.frame(values = values, sample = samples, alleles = alleles, na=na)
    df$sample <- factor(df$sample, levels = unique(df$sample),ordered = TRUE)
    df$alleles <- factor(df$alleles, levels = unique(df$alleles),ordered = TRUE)

	#if phase option add phase information to df
	if(top.allele.criteria=="phase"){
		df$phase <- phase
		df$phase <- factor(df$phase, levels = unique(df$phase),ordered = TRUE)
	}
	df

})

#' defaultPhase
#' 
#' used to populate the phase slot in an ASEset object
#' 
#' will set everything to 0
#'
#' @name defaultPhase
#' @rdname defaultPhase
#' @aliases defaultPhase,numeric-method
#' @docType methods
#' @param i number of rows
#' @param j number of columns 
#' @param ... arguments to forward to internal functions
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords phase
#' @examples
#' 
#' 
#' i <- 5
#' j <- 10
#' defaultPhase(i,j)
#'
NULL

#' @rdname defaultPhase
#' @export
setGeneric("defaultPhase", function(i, ... ){
    standardGeneric("defaultPhase")
})


#' @rdname defaultPhase
#' @export
setMethod("defaultPhase", signature("numeric"),
		function(i, j, ...
	){
	#x is rows
	#y is columns

    p1 <- matrix(rep(0, i*j), nrow=i, ncol=j)
    p2 <- matrix(rep(0, i*j), nrow=i, ncol=j)
    matrix(paste(p1,rep("/", i*j), p2, sep=""),i, j)

	
})


