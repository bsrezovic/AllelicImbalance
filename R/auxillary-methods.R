#'@include ASEset-class.R
NULL

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
#' top.fraction.criteria can take three options, maxcount, ref and phase. The top
#' allele will be every second row in the data frame, with start from row 2. 
#' The maxcount argument will put the allele with most reads on top of the 
#' bivariate fraction. Similarly the ref argument will put always the reference
#' allele on top. The phase arguments puts the maternal phase always on top.
#' The top.fraction.criteria for the ref or phase arguments requires that both ref
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
#' @param top.fraction.criteria 'maxcount', 'ref' or 'phase'
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
setGeneric("fractionPlotDf", function(x, snp,  strand="*", top.fraction.criteria="maxcount", ... 
	){
    standardGeneric("fractionPlotDf")
})

#' @rdname fractionPlotDf
#' @export
setMethod("fractionPlotDf", signature(x = "ASEset"),
		function(x, snp,  strand="*", top.fraction.criteria="maxcount", ...
	){

	##top.fraction.criteria
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

	afraction <- fraction(x[snprow], strand = strand, top.fraction.criteria=top.fraction.criteria)

	values <- as.vector(t(matrix(as.numeric(afraction), ncol=2, nrow=ncol(x))))
	#values[seq(1,ncol(x)*2,by=2)+1] <- 1 - values[seq(1,ncol(x)*2,by=2)+1]

	#The first value will be the bottom allele (second allele the top allele)
	values[seq(1,ncol(x)*2,by=2)] <- 1 - values[seq(1,ncol(x)*2,by=2)]
	samples <- as.vector(t(matrix(rownames(afraction), ncol=2, nrow=ncol(x))))

	if(top.fraction.criteria=="phase"){

		if(!sum(c("ref", "alt") %in% names(mcols(x))) == 2 ){
			stop("ref and alt has to be set in mcols to use phase option")
		}

		mat <- phase(x[snprow],return.class="array")[,,1]
		mat2 <- matrix(c(rep(mcols(x[snprow])[,"ref"],ncol(x)),rep(mcols(x[snprow])[,"alt"], ncol(x))), ncol=2, nrow=ncol(x))
		mat2[mat==0,1] <- mcols(x[snprow])[,"alt"]
		mat2[mat==0,2] <- mcols(x[snprow])[,"ref"]

		alleles <- as.vector(t(mat2))
		phase <- rep(c("paternal","maternal"), ncol(x))

	}else if(top.fraction.criteria=="ref"){
		
		if(!sum(c("ref", "alt") %in% names(mcols(x))) == 2 ){
			stop("ref and alt has to be set in mcols to use phase option")
		}

		alleles <- as.vector(t(matrix(c(rep(mcols(x[snprow])[,"alt"],ncol(x)),rep(mcols(x[snprow])[,"ref"], ncol(x))), ncol=2, nrow=ncol(x))))
		
	}else if(top.fraction.criteria=="maxcount"){

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
	if(top.fraction.criteria=="phase"){
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
#' #add alternative allele information
#' mcols(a)[["alt"]] <- inferAltAllele(a)
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

		fr <- fraction(x, strand=strand, top.fraction.criteria="phase")

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

#' scanForHeterozygotes
#' 
#' Identifies the positions of SNPs found in BamGR reads.
#' 
#' This function scans all reads stored in a \code{GAlignmentsList} for
#' possible heterozygote positions. The user can balance the sensitivity of the
#' search by modifying the minimumReadsAtPos, maximumMajorAlleleFrequency and
#' minimumBiAllelicFrequency arguments.
#' 
#' @name ASEset-scanForHeterozygotes
#' @rdname ASEset-scanForHeterozygotes
#' @aliases ASEset-scanForHeterozygotes scanForHeterozygotes scanForHeterozygotes,ASEset-method
#' @docType methods
#' @param BamList A \code{GAlignmentsList object}
#' @param minimumReadsAtPos minimum number of reads required to call a SNP at a
#' given position
#' @param maximumMajorAlleleFrequency maximum frequency allowed for the most
#' common allele. Setting this parameter lower will minimise the SNP calls
#' resulting from technical read errors, at the cost of missing loci with
#' potential strong ASE
#' @param minimumMinorAlleleFrequency minimum frequency allowed for the second most
#' common allele. Setting this parameter higher will minimise the SNP calls
#' resulting from technical read errors, at the cost of missing loci with
#' potential strong ASE
#' @param minimumBiAllelicFrequency minimum frequency allowed for the first and
#' second most common allele. Setting a Lower value for this parameter will
#' minimise the identification of loci with three or more alleles in one
#' sample. This is useful if sequencing errors are suspected to be common.
#' @param verbose logical indicating if process information should be displayed
#' @param ... argument to pass on
#' @return \code{scanForHeterozygotes} returns a GRanges object with the SNPs
#' for the BamList object that was used as input.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{getAlleleCounts}} which is a
#' function that count the number of reads overlapping a site.  }
#' @keywords scan SNP heterozygote
#' @examples
#' 
#' data(reads)
#' s <- scanForHeterozygotes(reads,verbose=FALSE)
#' 
NULL

#' @rdname ASEset-scanForHeterozygotes
#' @export
setGeneric("scanForHeterozygotes", function(BamList, ...){
    standardGeneric("scanForHeterozygotes")
})

#' @rdname ASEset-scanForHeterozygotes
#' @export
setMethod("scanForHeterozygotes", signature(BamList = "GAlignmentsList"), 
		  function(BamList, minimumReadsAtPos = 20,
			maximumMajorAlleleFrequency = 0.90, minimumMinorAlleleFrequency = 0.1, 
			minimumBiAllelicFrequency = 0.9, verbose = TRUE, ...) {
    
    # if just one element of, make list (which is a convenient way of handling this
    # input type)
    if (class(BamList) == "GAlignments") {
        BamList <- GAlignmentsList(BamList)
    }
    if (class(BamList) == "GRanges") {
        BamList <- GRangesList(BamList)
    }
    
    if (class(minimumReadsAtPos) != "numeric") 
        stop(paste("minimumReadsAtPos must be of class numeric, not", class(minimumReadsAtPos)))
    if (length(minimumReadsAtPos) != 1) 
        stop(paste("minimumReadsAtPos must be of length 1, not", length(minimumReadsAtPos)))
    if (class(maximumMajorAlleleFrequency) != "numeric") 
        stop(paste("maximumMajorAlleleFrequency must be of class numeric, not", class(maximumMajorAlleleFrequency)))
    if (length(maximumMajorAlleleFrequency) != 1) 
        stop(paste("maximumMajorAlleleFrequency must be of length 1, not", length(maximumMajorAlleleFrequency)))
    if (maximumMajorAlleleFrequency < 0 | maximumMajorAlleleFrequency > 1) 
        stop("maximumMajorAlleleFrequency must be between 0 and 1")
    if (class(minimumBiAllelicFrequency) != "numeric") 
        stop(paste("minimumBiAllelicFrequency must be of class numeric, not", class(minimumBiAllelicFrequency)))
    if (length(minimumBiAllelicFrequency) != 1) 
        stop(paste("minimumBiAllelicFrequency must be of length 1, not", length(minimumBiAllelicFrequency)))
    if (minimumBiAllelicFrequency < 0 | maximumMajorAlleleFrequency > 1) 
        stop("minimumBiAllelicFrequency must be between 0 and 1")
    if (class(verbose) != "logical") 
        stop(paste("verbose must be of class logical, not", class(verbose)))
    if (length(verbose) != 1) 
        stop(paste("verbose must be of length 1, not", length(verbose)))
    
    # checking that BamList format is ok (has to be done quite carefully for the
    # list-class gappedAlignments
    if (class(BamList) == "GRangesList") 
        stop("The use of GRangesList is not recommended. Use BamImpGAList or BamImpGAPList")
    if (class(BamList) != "GAlignmentsList") 
        stop("The class of the BamList has to be a GAlignmentsList")
    
    # checks specific for GAlignments
    if (unique(unlist(lapply(BamList, class))) == "GAlignments") {
        if (!all(unlist(lapply(BamList, function(x) {
            all(c("cigar", "qwidth") %in% colnames(mcols(x)))
        })))) {
            stop("BamList given as GAlignmentsLists of GAlignments objects must contain mcols named qwidth and cigar")
        }
    }
    
    
    RangedData <- GRanges()
    chromosomeLevels <- unique(unlist(lapply(BamList, function(x) {
        levels(droplevels(runValue(seqnames(x))))
    })))
    
	#create pos to extract
	x <- BamList
	x <- GAlignmentsList(lapply(x,function(y){
			strand(y)<- "*"
			y	
	}))
	x <- reduce(granges(unlist(x)))
	pos <- unlist(mapply(width(x), start(x), FUN = function(y, z){z:(z+y-1)}))
	chrnames <- unlist(mapply(as.character(seqnames(x)),width(x), FUN=rep))
	#strand.vec <- unlist(mapply(as.character(strand(x)),width(x), FUN=rep))
	my_IGPOI <- GRanges(seqnames=chrnames, ranges=IRanges(start=pos, width=1), 
						strand="*")

	#empty array that handles only four nucleotides + one del columns
	dimnames = list(paste("pos",1:length(my_IGPOI), sep=""), names(BamList), c("A", "C", "G", "T"))
	ar <- array(NA, c(length(my_IGPOI), length(BamList), 4),
				 dimnames = dimnames)  

			
	for (i in 1:length(BamList)) {
		if (verbose) 
			cat(paste("Investigating sample", i), "out of", length(BamList), 
			  "\n")
		
		# extract samples
		gal <- BamList[[i]]
		
		if (!(length(gal) == 0)) {
                
			seqlevels(gal) <- seqlevels(my_IGPOI) 
			qseq <- mcols(gal)$seq

			nuclpiles <- pileLettersAt(qseq, seqnames(gal), start(gal), cigar(gal),
																 my_IGPOI)
			# fill array
			nstr <- strsplit(as.character(nuclpiles), "")
			for (k in 1:length(my_IGPOI)) {
				ar[k, i, ] <- c(sum(nstr[[k]] %in% "A"), sum(nstr[[k]] %in% "C"), sum(nstr[[k]] %in% 
					"G"), sum(nstr[[k]] %in% "T"))  
			}
		}
	}	

	#mat <- apply(ar, c(1,3),sum, na.rm=TRUE)
	allele.count.tot <- apply(ar, c(1,2), sum)
	tf <- allele.count.tot  < minimumReadsAtPos
	allele.count.tot[tf] <- NA


	#make frequency array
	ar2 <- ar * array(as.vector(1/allele.count.tot),dim=dim(ar))
	ar2[is.na(ar2)] <- 0

	#rank alleles
	rank <- t(apply(apply(ar,c(1,3),sum, na.rm=TRUE),
					1, function(x){rank(x,ties.method="first")}))

	ar.rank1 <- array(rank==4, dim=c(nrow(ar2), dim(ar2)[3], ncol=ncol(ar2)))
	ar.rank2 <- array(rank==3, dim=c(nrow(ar2), dim(ar2)[3], ncol=ncol(ar2)))

	#rearrange to be able to subset
	ar.rank1 <- aperm(ar.rank1, c(1,3,2))
	ar.rank2 <- aperm(ar.rank2, c(1,3,2))

	#rearrange to be able to transform back
	ar3 <- aperm(ar2, c(3,2,1))
	ar.rank1b <- aperm(ar.rank1, c(3,2,1))
	ar.rank2b <- aperm(ar.rank2, c(3,2,1))

	#subset ref allele frequencies
	maj.al.fr <- aperm(array(ar3[ar.rank1b],dim=dim(ar2)[2:1], 
					   dimnames=dimnames(ar2)[2:1]),c(2,1))
					   
	min.al.fr <- aperm(array(ar3[ar.rank2b],dim=dim(ar2)[2:1], 
					   dimnames=dimnames(ar2)[2:1]),c(2,1))

	#check conditions
	tf1 <- apply(!min.al.fr <= minimumMinorAlleleFrequency,1,any)
	tf2 <- apply(!maj.al.fr >= maximumMajorAlleleFrequency,1,any)
	tf3 <- apply(((min.al.fr + maj.al.fr) >= minimumBiAllelicFrequency),1,any)

	toBeReturned <- my_IGPOI[tf1&tf2&tf3]
    
    #add an SNP name based on position
    if (!(length(toBeReturned) == 0)) {
        names(toBeReturned) <- paste("chr", seqnames(toBeReturned), "_", start(toBeReturned), 
            sep = "")
        
    }

    return(toBeReturned)
})

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
#' @aliases import-bam impBamGAL impBamGAL,character-method
#' @docType methods
#' @param UserDir The relative or full path of folder containing bam files.
#' @param searchArea A \code{GenomicRanges object} that contains the regions of
#' interest
#' @param XStag Setting \code{XStag=TRUE} stores the strand specific
#' information in the mcols slot 'XS'
#' @param verbose Setting \code{verbose=TRUE} gives details of procedure during
#' function run.
#' @param ... arguments to pass on
#' @return \code{impBamGRL} returns a GRangesList object containing the RNA-seq
#' reads in the region defined by the \code{searchArea} argument.
#' \code{impBamGAL} returns a list with GAlignments objects containing the
#' RNA-seq reads in the region defined by the \code{searchArea} argument.
#' \code{funImpBamGAPL} returns a list with GappedAlignmentPairs object
#' containing the RNA-seq reads in the region defined by the \code{searchArea}
#' argument.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords bam import
#' @examples
#' 
#' #Declare searchArea
#' searchArea <- GRanges(seqnames=c('17'), ranges=IRanges(79478301,79478361))
#' 
#' #Relative or full path  
#' pathToFiles <- system.file('extdata/ERP000101_subset', package='AllelicImbalance')
#' 
#' reads <- impBamGAL(pathToFiles,searchArea,verbose=FALSE)
#' 
NULL

#' @rdname import-bam
#' @export
setGeneric("impBamGAL", function(UserDir, ...){
    standardGeneric("impBamGAL")
})

#' @rdname import-bam
#' @export
setMethod("impBamGAL", signature(UserDir = "character"), 
	function(UserDir, searchArea, XStag = FALSE, verbose = TRUE, ...) {
    # Set parameters
    which <- searchArea  #A GRanges, RangesList, RangedData, or missing object, from which a IRangesList instance will be constructed.
    what <- scanBamWhat()  #A character vector naming the fields to return. scanBamWhat() returns a vector of available fields. Fields are described on the scanBam help page.
    flag <- scanBamFlag(isUnmappedQuery = FALSE)
    
    if (XStag) {
        param <- ScanBamParam(flag = flag, which = which, what = what, tag = "XS")  #store ScanBamParam in param.
    } else {
        param <- ScanBamParam(flag = flag, which = which, what = what)  #store ScanBamParam in param.\t\t
    }
    # Point to correct directory and create a BamFileList object
    bamDir <- normalizePath(UserDir)  #Point to the directory containing your Bam files and its respective bam.bai files.
    allFiles <- list.files(bamDir, full.names = TRUE)  #list files in a folder.
    bamFiles <- allFiles[grep(".bam$", allFiles)]  #list only the files ending in .bam .
    if (length(bamFiles) == 0) 
        stop(paste("No bam files found in", bamDir))
    if (!all(file.exists(paste(bamFiles, ".bai", sep = "")))) {
        if (verbose) 
            cat(paste("The bam files in UserDir are required to also have .bam.bai index files. Trying to run indexBam function on each"), 
                "\n")
        indexBam(bamFiles)
        if (!all(file.exists(paste(bamFiles, ".bai", sep = "")))) {
            stop("The bam files in UserDir are required to also have .bam.bai index files.")
        } else {
            if (verbose) 
                cat(paste("Succesfully indexed all bamFiles in UserDir", UserDir), 
                  "\n")
        }
    }
    bamFilesList <- BamFileList(bamFiles)  #store all the .bam paths in a BamFile.
    
    # check that sequences in searchArea are actually found in the bam files
    header <- scanBamHeader(bamFiles)
    checkSeqNameExists <- function(bamHeader, requestedSeqNames) {
        as.character(requestedSeqNames) %in% names(bamHeader[["targets"]])
    }
    if (!all(unlist(lapply(header, checkSeqNameExists, seqnames(searchArea))))) {
        # not all searchArea requested seq-names found in bam files. Create nice error
        # report and stop
        seqNotFoundErrors <- lapply(header, checkSeqNameExists, seqnames(searchArea))
        seqNotFounds <- vector()
        for (sampleName in names(seqNotFoundErrors)) {
            seqNotFounds <- c(seqNotFounds, as.character(seqnames(searchArea)[!seqNotFoundErrors[[sampleName]]]))
        }
        stop(paste("The following seq name(s) not found in the bam files:", paste(sort(unique(seqNotFounds)), 
            collapse = ", ")))
    }
    
    # Loop through, open scanBam, store in GRList and then close each object in the
    # BamFileList object.
    BamGAL <- list()
    i <- 1
    for (bamName in names(bamFilesList)) {
        # Description
        bf <- bamFilesList[[bamName]]
        open(bf)
        if (verbose) 
            cat(paste("Reading bam file", i, "with filename", basename(bamName)), 
                "\n")  #Print information to the user
        GappedAlign <- readGAlignments(bf, param = param)
        
        BamGAL[[basename(bamName)]] <- GappedAlign
        
        if (verbose) 
            cat(paste("stored", basename(bamName), "in BamGAL"), "\n")
        gc()
        close(bf)
        i <- i + 1
    }
    BamGAL <- GAlignmentsList(BamGAL)
    
    return(BamGAL)
})



#' Import Bcf Selection
#' 
#' Imports a selection of a bcf file or files specified by a GenomicRanges
#' object as search area.
#' 
#' A wrapper to import bcf files into R in the form of GenomicRanges objects.
#' 
#' @name import-bcf
#' @rdname import-bcf
#' @aliases import-bcf impBcfGRL impBcfGR impBcfGRL,character-method impBcfGR,character-method
#' @docType methods
#' @param UserDir The relative or full path of folder containing bam files.
#' @param searchArea A \code{GenomicRanges} object that contains the regions of
#' interest
#' @param verbose Setting \code{verbose=TRUE} gives details of the procedure
#' during function run.
#' @param ... parameters to pass on
#' @return \code{BcfImpGRList} returns a GRangesList object.  \code{BcfImpGR}
#' returns one GRanges object of all unique entries from one or more bcf files.
#' @note Make sure there is a complementary index file \code{*.bcf.bci} for
#' each bcf file in \code{UserDir}. If there is not, then the functions
#' \code{impBcfGRL} and \code{impBcfGR} will try to create them.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The impBamGRL for importing bam files
#' \item The \code{\link{getAlleleCounts}} for how to get allele(SNP) counts
#' \item The \code{\link{scanForHeterozygotes}} for how to find possible
#' heterozygote positions }
#' @keywords bcf import
#' @examples
#' 
#' #Declare searchArea
#' searchArea <- GRanges(seqnames=c('17'), ranges=IRanges(79478301,79478361))
#' 
#' #Relative or full path  
#' pathToFiles <- system.file('extdata/ERP000101_subset', package='AllelicImbalance')
#' 
#' #import
#' reads <- impBcfGRL(pathToFiles, searchArea, verbose=FALSE)
#' 
NULL

#' @rdname import-bcf
#' @export
setGeneric("impBcfGRL", function(UserDir, ... 
	){
    standardGeneric("impBcfGRL")
})

#' @rdname import-bcf
#' @export
setMethod("impBcfGRL", signature(UserDir = "character"),
	function(UserDir, searchArea = NULL, verbose = TRUE, ...) {
    
    # Set parameters
    if (is.null(searchArea)) {
        param <- ScanBcfParam()
    } else {
        param <- ScanBcfParam(which = searchArea)
    }
    # Point to correct directory and create a BcfFileList object
    bcfDir <- normalizePath(UserDir)  #Point to the directory containing your Bam files and its respective bam.bai files.
    allFiles <- list.files(bcfDir, full.names = TRUE)  #list files in a folder.
    bcfFiles <- allFiles[grep(".bcf$", allFiles)]  #list only the files ending in .bam .
    if (length(bcfFiles) == 0) 
        stop(paste("No bcf files were found in", UserDir))
    
    # bcfFilesList <- BcfFileList(bcfFiles) #store all the .bam paths in a BamFile.
    if (!all(file.exists(paste(bcfFiles, ".bci", sep = "")))) {
        if (verbose) 
            cat("Did not find bci files for all bcf files. Trying the indexBcf function obtain these", 
                "\n")
        for (bcfFile in bcfFiles) {
            indexBcf(bcfFile)
        }
        if (!all(file.exists(paste(bcfFiles, ".bci", sep = "")))) {
            stop("The bcf files in UserDir are required to also have .bcf.bci index files. Run the indexBcf function in package Rsamtools on each bam file.")
        }
    }
    
    # Loop through, open scanBam, store in GRList and then close each object in the
    # BamFileList object.
    BcfGRL <- GRangesList()
    for (i in 1:length(bcfFiles)) {
        
        bcf <- suppressWarnings(scanBcf(file = bcfFiles[i], param = param))
        
        # need to protect against empty bcf files
        if (length(bcf[["POS"]]) == 0) {
            GRangeBcf <- GRanges(seqnames = vector(), ranges = IRanges(start = vector(), 
                width = vector()), ref = vector(), alt = vector(), qual = vector())
            bcfName <- bcfFiles[i]
            BcfGRL[[basename(bcfName)]] <- GRangeBcf
            
        } else {
            # if they are not empty we just proceed as usual
            ranges <- IRanges(start = bcf[["POS"]], width = 1L)
            GRangeBcf <- GRanges(seqnames = as.character(bcf[["CHROM"]]), ranges = ranges, 
                ref = bcf[["REF"]], alt = bcf[["ALT"]], qual = bcf[["QUAL"]])
            # Store GRangeBam in BamGRL (which is the GRange List object)
            bcfName <- bcfFiles[i]
            BcfGRL[[basename(bcfName)]] <- GRangeBcf
            if (verbose) 
                cat(paste("stored", basename(bcfName), "in BcfGRL"), "\n")
            gc()
        }
    }
    return(BcfGRL)
})

#' @rdname import-bcf
#' @export
setGeneric("impBcfGR", function(UserDir, ... 
	){
    standardGeneric("impBcfGR")
})

#' @rdname import-bcf
#' @export
setMethod("impBcfGR", signature(UserDir = "character"),
	function(UserDir, searchArea = NULL, verbose = TRUE, ...) {
		BcfGRList <- impBcfGRL(UserDir, searchArea, verbose)
		BcfGR <- do.call(c, unname(as.list(BcfGRList)))
		BcfGR <- unique(BcfGR)
		names(BcfGR) <- paste("chr", seqnames(BcfGR), "_", start(BcfGR), sep = "")
		return(BcfGR)
})

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
#' corresponding to bi-allelic SNPs. In the \code{strand} argument, specifying
#' '*' is the same as retrieving the sum count of '+' and '-' reads
#' (and unknown strand reads in case these are found in the bam file). '*' is
#' the default behaviour and can be used when the RNA-seq experiments strand
#' information is not available.
#' 
#' @name getAlleleCounts
#' @rdname getAlleleCounts
#' @aliases BamList getAlleleCounts getAlleleCounts,GAlignmentsList-method
#' @docType methods
#' @param BamList A \code{GAlignmentsList object} or \code{GRangesList object}
#' containing data imported from a bam file
#' @param GRvariants A \code{GRanges object} that contains positions of SNPs to
#' retrieve
#' @param strand A length 1 \code{character} with value  '+',
#' '-', or '*'.  This argument determines if \code{getAlleleCounts} will
#' retrieve counts from all reads, or only from reads marked as '+', '-' or '*'
#' (unknown), respectively.
#' @param return.class 'list' or 'array'
#' @param verbose Setting \code{verbose=TRUE} makes function more talkative
#' @param ... parameters to pass on
#' @return \code{getAlleleCounts} returns a list of several data.frame objects,
#' each storing the count data for one SNP.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{scanForHeterozygotes}} which is a
#' function to find possible heterozygote sites in a
#' \code{\link[GenomicAlignments]{GAlignmentsList}} object }
#' @keywords SNP count
#' @examples
#' 
#' #load example data
#' data(reads)
#' data(GRvariants)
#' 
#' 
#' #get counts at the three positions specified in GRvariants
#' alleleCount <- getAlleleCounts(BamList=reads,GRvariants,
#' strand='*')
#' 
#' #if the reads had contained stranded data, these two calls would 
#' #have given the correct input objects for getAlleleCounts
#' alleleCountPlus <- getAlleleCounts(BamList=reads,GRvariants,
#' strand='+')
#' alleleCountMinus <- getAlleleCounts(BamList=reads,GRvariants,
#' strand='-')
#' 
#' 
NULL

#' @rdname getAlleleCounts
#' @export
setGeneric("getAlleleCounts", function(BamList, ... 
	){
    standardGeneric("getAlleleCounts")
})

#' @rdname getAlleleCounts
#' @export
setMethod("getAlleleCounts", signature(BamList = "GAlignmentsList"),
function(BamList, GRvariants, strand = "*",
						return.class = "list", verbose = TRUE, ...) { 
    
    
    if (!class(BamList) %in% c("GAlignments", "GAlignmentsList")) {
        stop("BamList has to be of class GAlignments or GAlignmnetsList\n")
    }
    # if just one element of, make list (which is a convenient way of
	# handling this input type)
    # 
    if (class(BamList) == "GAlignments") {
        BamList <- GAlignmentsList(BamList)
    }
    
    # check for strand name
    if (!class(strand) == "character") {
        stop("strand has to be of class character")
    }
    if (!length(strand) == 1) {
        stop("strand has to be of length 1")
    }
    if (!sum(strand %in% c("+", "-", "*")) > 0) {
        stop("strand parameter has to be either '+', '-' or '*' ")
    }
    
    # if the user sent in the GRangesList for GRvariants,
	# take out only the unique entries.
    # 
    if (class(GRvariants) == "GRangesList") {
        GRvariants <- unique(unlist(GRvariants, use.names = FALSE)) 
    }
   
	#if BamList is not list, make it a list
	if(class(BamList)=="GAlignments"){
		BamList <- GAlignmentsList(BamList)
	}

	#Drop seqlevels in BamList that are not in GRvariants
	#seqlevels(BamList,force=TRUE) <- seqlevels(GRvariants)
	seqinfo(GRvariants) <- merge(seqinfo(GRvariants), seqinfo(BamList))
	seqlevels(GRvariants) <- seqlevelsInUse(GRvariants)


    # check that seqlevels are the same
   # if (!identical(seqlevels(BamList), seqlevels(GRvariants))) {
   #     stop("!identical(seqlevels(BamList), seqlevels(GRvariants))\n")
   # }
    
    # checking that GRvariants is ok
    if (class(GRvariants) != "GRanges") 
        stop(paste("GRvariants must be of class GRanges, not",
				   class(GRvariants)))
    if (length(GRvariants) == 0) 
        stop("GRvariants was given as an empty GRanges object.",
			 " There can be no Snps retrieved by getAlleleCount then")
    if (any(width(GRvariants) != 1)) 
        stop("GRvariants can contain only entries of width=1,",
			 " corresponding to SNPs.")
    
    # checking that verbose is ok
    if (class(verbose) != "logical") 
        stop(paste("verbose must be of class logical, not", class(verbose)))
    if (length(verbose) != 1) 
        stop(paste("verbose must be of length 1, not", length(verbose)))
    
    # make row-names
    if (sum(grepl("chr", seqnames(GRvariants))) > 0) {
        snpNames <- paste(seqnames(GRvariants),
						  "_", start(GRvariants), sep = "")
    } else {
        snpNames <- paste("chr", seqnames(GRvariants),
						  "_", start(GRvariants), sep = "")
    }
    
	# needs name, need a more general solution here
	if(length(names(BamList)) == 0){
		warning("no set names for list, new names will be sample1,2,3,etc")
		names(BamList) <- paste("sample",1:length(BamList),sep="")
	}

	#empty array that handles only four nucleotides + one del columns
    dimnames = list(snpNames, names(BamList), c("A", "C", "G", "T"))
    ar1 <- array(NA, c(length(GRvariants), length(BamList), 4),
				 dimnames = dimnames)  
    
    # use strand choice to only get reads from that strand
    if (!strand == "*") {
        BamList <- GAlignmentsList(mapply(function(x, y) {
            x[y]
        }, BamList, strand(BamList) == strand))
    }
    
    for (j in 1:length(names(BamList))) {
        sample <- names(BamList)[j]
        if (verbose) 
            cat("sample ", sample, "\n")
        
        gal <- BamList[[j]]
		my_IGPOI <- GRvariants
		seqlevels(gal) <- seqlevels(my_IGPOI) 
		qseq <- mcols(gal)$seq
        
#        nuclpiles <- pileLettersAt(mcols(gal)[, "seq"], seqnames(gal), start(gal), 
#            cigar(gal), GRvariants)
#
		nuclpiles <- pileLettersAt(qseq, seqnames(gal), start(gal), cigar(gal),
							                                 my_IGPOI)
        
        # fill array
        nstr <- strsplit(as.character(nuclpiles), "")
        for (k in 1:length(GRvariants)) {
            ar1[k, j, ] <- c(sum(nstr[[k]] %in% "A"), sum(nstr[[k]] %in% "C"), sum(nstr[[k]] %in% 
                "G"), sum(nstr[[k]] %in% "T"))  
        }
        
    }
    
    # check return.type argument
    if (return.class == "list") {
        alleleCountList <- list()
        for (i in 1:nrow(ar1)) {
            mat <- ar1[i, , ]
			if(class(mat)=="integer"){
				mat <- t(as.matrix(mat))
			}
            if (class(mat) == "numeric") {
                mat <- t(mat)
                colnames(mat) <- dimnames[[3]]
            } else {
                colnames(mat) <- dimnames[[3]]
            }
            rownames(mat) <- dimnames[[2]]
            alleleCountList[[i]] <- mat
        }
        names(alleleCountList) <- dimnames[[1]]
        alleleCountList
    } else if (return.class == "array") {
        ar1
    } else {
        cat("return.type unknown\n Nothing will be returned from function!")
    }
})

#' Map Bias
#' 
#' an allele frequency array
#' 
#' This function will assume there is no bias that comes from the mapping of
#' reads, and therefore create a matrix with expected frequency of 0.5 for each
#' allele.
#' 
#' @name getDefaultMapBiasExpMean
#' @rdname getDefaultMapBiasExpMean
#' @aliases getDefaultMapBiasExpMean getDefaultMapBiasExpMean3D 
#' getDefaultMapBiasExpMean,ANY-method getDefaultMapBiasExpMean3D,ANY-method
#' @aliases getDefaultMapBiasExpMean getDefaultMapBiasExpMean3D
#' @docType methods
#' @param alleleCountList A \code{GRangesList object} containing read
#' information
#' @param ... parameters to pass on
#' @return \code{getDefaultMapBiasExpMean} returns a matrix with a default
#' expected mean of 0.5 for every element.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords mapping bias
#' @examples
#' 
#' #load example data
#' data(ASEset)
#' #access SnpAfList
#' alleleCountList <- alleleCounts(ASEset)
#' #get default map bias exp mean
#' matExpMean <- getDefaultMapBiasExpMean(alleleCountList)
#' 
#' 
NULL

#' @rdname getDefaultMapBiasExpMean
#' @export
setGeneric("getDefaultMapBiasExpMean", function(alleleCountList, ... 
	){
    standardGeneric("getDefaultMapBiasExpMean")
})

#' @rdname getDefaultMapBiasExpMean
#' @export
setGeneric("getDefaultMapBiasExpMean3D", function(alleleCountList, ... 
	){
    standardGeneric("getDefaultMapBiasExpMean3D")
})

#' @rdname getDefaultMapBiasExpMean
#' @export
setMethod("getDefaultMapBiasExpMean", signature(alleleCountList = "list"),
 function(alleleCountList) {
    
    l <- lapply(alleleCountList, function(x) {
        ap <- apply(x, 2, sum)
        char <- names(sort(ap, decreasing = TRUE))[1:2]
        
        v <- rep(0, length(colnames(x)))
        v[colnames(x) %in% char] <- 0.5
        v
    })
    
    MapBiasExpMean <- matrix(unlist(l), byrow = TRUE, nrow = length(alleleCountList), 
        ncol = 4, dimnames = list(c(names(alleleCountList)), colnames(alleleCountList[[1]])))  # alleleCountList[[1]] assumes that in each list the colnames are the same.
    MapBiasExpMean
})


#' @rdname getDefaultMapBiasExpMean
#' @export
setMethod("getDefaultMapBiasExpMean3D", signature(alleleCountList = "ANY"),
function(alleleCountList) {
   
	if(class(alleleCountList)=="list"){
		MapBiasExpMean <- getDefaultMapBiasExpMean(alleleCountList)
		# make 3D array
		MapBiasExpMean3D <- array(NA, c(length(alleleCountList),
										length(unlist(unique(lapply(alleleCountList, 
			rownames)))), 4))  #empty array
		for (i in 1:length(unlist(unique(lapply(alleleCountList, rownames))))) {
			MapBiasExpMean3D[, i, ] <- MapBiasExpMean
		}
	}
	if(class(alleleCountList)=="array"){
		# make 3D array
		MapBiasExpMean3D <- alleleCountList
		mapbiasmat <- t(apply(apply(alleleCountList,c(1,3),sum),
					1, function(x){rank(x,ties.method="first")}))						
		mapbiasmat[mapbiasmat==1] <- 0.5
		mapbiasmat[mapbiasmat==2] <- 0.5
		mapbiasmat[mapbiasmat==3] <- 0
		mapbiasmat[mapbiasmat==4] <- 0

		MapBiasExpMean3D <- array(NA, dim=c(nrow(mapbiasmat),dim(alleleCountList)[2],4))
		for (i in 1:dim(alleleCountList)[2]) {
			MapBiasExpMean3D[, i, ] <- mapbiasmat
		}

	}
    MapBiasExpMean3D
})

