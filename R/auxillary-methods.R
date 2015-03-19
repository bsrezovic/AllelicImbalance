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
		phased <- as.integer(upsplit[seq.int(from=2, to=length(upsplit), by=3)]=="|")
	
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


#' regionSummary
#' 
#' Gives a summary of AI-consistency for a transcript
#'
#' From a given set of e.g. transcripts exon ranges the function will return
#' a summary for the sum of all exons. Phase information is required.
#'
#' @name regionSummary
#' @rdname regionSummary
#' @aliases regionSummary,numeric-method
#' @docType methods
#' @param x ASEset object
#' @param strand can be "+", "-" or "*"
#' @param gr GenomicRanges object to summmarize over
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
#' # in this example every snp is on its own exon
#' txGR <- granges(a)
#' t <- regionSummary(a, txGR)
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
		function(x, gr, strand="*", ...
	){

		#needs alternative allele
		if(!"alt" %in% colnames(mcols(x))){
			stop("function needs mcols(x)[['alt']] to be set")				
		}

		#make overlap and subset based on gr
		hits <- findOverlaps(x,gr)
		x <- x[queryHits(hits),]

		fr <- fraction(x, strand=strand, top.allele.criteria="phase")

		#need information of which are heterozygotes and homozygotes
		if(is.null(genotype(x))){
			genotype(x) <- inferGenotypes(x, return.allele.allowed="bi")
		}
		fr.het.filt <- hetFilt(x)
		fr.f <- fr
		fr.f[!t(fr.het.filt)] <- NaN

		maternalAllele <- function(x){
		
			mat <- phase(x,return.class="array")[,,1]
			ref <- mcols(x)[["ref"]]
			alt <- mcols(x)[["alt"]]
		
			apply(t(mat),1,function(y){
				
				vec <- rep(NA,length(y))
				if(any(y == 1)){
					vec[y == 1] <- ref[y == 1]
				}
				if(any(y == 0)){
					vec[y == 0] <- alt[y == 0]
				}
				vec
			})
			
		}

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

		tf <- fr.f > t(mbias.values)
		tf2 <- fr.f < t(mbias.values)
		fr.up.filt <- tf
		fr.down.filt <- tf2

		pv <- binom.test(x,strand)

		p.up.filt <- pv
		p.up.filt[!tf] <- NA
		p.down.filt <- pv
		p.down.filt[!tf2] <- NA

		ai.down <- apply(tf2,1,sum, na.rm=TRUE)
		ai.up <- apply(tf,1,sum, na.rm=TRUE)

		fr.d <- abs(fr.f - t(mbias.values))

		hets <- apply(t(fr.het.filt),1,sum)
		homs <- apply(t(!fr.het.filt),1,sum)

		mean.fr <- apply(fr.f, 1, mean, na.rm=TRUE)
		sd.fr <- apply(fr.f, 1, sd, na.rm=TRUE)
		mean.delta <- apply(fr.d, 1, mean, na.rm=TRUE)
		sd.delta <- apply(fr.d, 1, sd, na.rm=TRUE)

		#dir.up <- apply(fr.f,1,function(x){sum(mean(x, na.rm=TRUE)>0.5)})
		#dir.down <- apply(fr.f,1,function(x){sum(mean(x, na.rm=TRUE)<0.5)})

		#return data frame
		data.frame(
				het=hets,
				hom=homs,
				mean.fr=mean.fr,
				sd.fr=sd.fr,
				mean.delta=mean.delta,
				sd.delta=sd.delta,
				ai.up=ai.up,
				ai.down=ai.down
				)
})



