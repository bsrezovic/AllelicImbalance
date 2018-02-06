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
#' @param dimnames list with dimnames
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
setGeneric("phaseMatrix2Array", function(x, ...)
    standardGeneric("phaseMatrix2Array")
)

#' @rdname phaseMatrix2Array
#' @export
setMethod("phaseMatrix2Array", signature(x = "matrix"),
		function(x, dimnames=NULL, ...
	){

		psplit <- strsplit(x, split="")
		upsplit <- unlist(psplit)
		mat <- as.integer(upsplit[seq(1, length(upsplit), by=3)])
		pat <- as.integer(upsplit[seq(3, length(upsplit), by=3)])
		phased <- as.integer(upsplit[seq.int(from=2, to=length(upsplit), by=3)]=="|")

		array(c(mat, pat, phased), dim=c(nrow(x), ncol(x), 3),
			  dimnames = dimnames)

})

#' phaseArray2phaseMatrix
#'
#' used to convert the phase from the visually friendly matrix to array.
#'
#' A more effectice way of store the phase data in the ASEset object
#'
#' @name phaseArray2phaseMatrix
#' @rdname phaseArray2phaseMatrix
#' @aliases phaseArray2phaseMatrix,array-method
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
#' mat <- phaseArray2phaseMatrix(ar)
#'
NULL

#' @rdname phaseArray2phaseMatrix
#' @export
setGeneric("phaseArray2phaseMatrix", function(x, ...)
	{
		standardGeneric("phaseArray2phaseMatrix")
	}
)

#' @rdname phaseArray2phaseMatrix
#' @export
setMethod("phaseArray2phaseMatrix", signature(x = "array"),
	function(x, ...)
	{
		.mergePhaseArray2phaseMatrix(x)
	}
)

### -------------------------------------------------------------------------
### helpers for phaseArray2phaseMatrix
###

# if any of the maternal paternal alles are NA, NA is retrieved
# if phase is NA it will be set to /
.mergePhaseArray2phaseMatrix <- function(x, ...)
{
	pha <- matrix("/", ncol=ncol(x), nrow=nrow(x))
	pha[x[,,3]==1] <- "|"

	nas <- is.na(x[,,1]) | is.na(x[,,2])
	merg <- paste(x[,,1], pha, x[,,2], sep="")
	merg[nas] <- NA

	matrix(merg, nrow=nrow(x), ncol(x))
}


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
	)
    standardGeneric("fractionPlotDf")
)

#' @rdname fractionPlotDf
#' @export
setMethod("fractionPlotDf", signature(x = "ASEset"),
		function(x, snp,  strand="*", top.fraction.criteria="maxcount", ...
	){

	##top.fraction.criteria
	# maxcount
	# ref
	# phase

	if(is.integer(snp)){
		snprow <- snp
	}else if(is.character(snp)){
		snprow <- which(rownames(x) %in% snp  )
	}else{
		# If the class of snp is not approprite there will be an informative error message from
		# the base subset method.
		snprow <- snp
	}

	afraction <- t(fraction(x[snprow], strand = strand, top.fraction.criteria=top.fraction.criteria))

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
    if (is(BamList, "GAlignments")) {
        BamList <- GAlignmentsList(BamList)
    }
    if (is(BamList, "GRanges")) {
        BamList <- GRangesList(BamList)
    }

    if (!is.numeric(minimumReadsAtPos))
        stop(paste("minimumReadsAtPos must be a numeric vector, not a", class(minimumReadsAtPos)))
    if (length(minimumReadsAtPos) != 1)
        stop(paste("minimumReadsAtPos must be of length 1, not", length(minimumReadsAtPos)))
    if (!is.numeric(maximumMajorAlleleFrequency))
        stop(paste("maximumMajorAlleleFrequency must be a numeric vector, not a", class(maximumMajorAlleleFrequency)))
    if (length(maximumMajorAlleleFrequency) != 1)
        stop(paste("maximumMajorAlleleFrequency must be of length 1, not", length(maximumMajorAlleleFrequency)))
    if (maximumMajorAlleleFrequency < 0 | maximumMajorAlleleFrequency > 1)
        stop("maximumMajorAlleleFrequency must be between 0 and 1")
    if (!is.numeric(minimumBiAllelicFrequency))
        stop(paste("minimumBiAllelicFrequency must be a numeric vector, not a", class(minimumBiAllelicFrequency)))
    if (length(minimumBiAllelicFrequency) != 1)
        stop(paste("minimumBiAllelicFrequency must be of length 1, not", length(minimumBiAllelicFrequency)))
    if (minimumBiAllelicFrequency < 0 | maximumMajorAlleleFrequency > 1)
        stop("minimumBiAllelicFrequency must be between 0 and 1")
    if (!is.logical(verbose))
        stop(paste("verbose must be of a logical, not a", class(verbose)))
    if (length(verbose) != 1)
        stop(paste("verbose must be of length 1, not", length(verbose)))

    # checking that BamList format is ok (has to be done quite carefully for the
    # list-class gappedAlignments
    if (is(BamList, "GRangesList"))
        stop("The use of GRangesList is not recommended. Use BamImpGAList or BamImpGAPList")
    if (!is(BamList, "GAlignmentsList"))
        stop("BamList has to be a GAlignmentsList object")

    # checks specific for GAlignments
    if (all(vapply(BamList, is, logical(1), "GAlignments"))) {
        if (!all(unlist(lapply(BamList, function(x) {
            all(c("cigar", "qwidth") %in% colnames(mcols(x)))
        })))) {
            stop("BamList given as GAlignmentsLists of GAlignments objects must contain mcols named qwidth and cigar")
        }
    }


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
#' Imports a specified genomic region from a bam file using a GRanges
#' object as search area.
#'
#' If the sequence data is strand-specific you may want to set XStag=TRUE. The
#' strand specific information has then to be stored in the meta columns with
#' column name 'XS'. If the aligner did not set the XS-tag and the data is strand-
#' specific it is still be possible to infer the strand from the bit flags after importing
#' the reads to R. Depending on the strand-specific protocol different combinations of the
#' flags will have to be used. For illumina fr-secondstrand, 83 and 163 are minus strand
#' reads and 99 and 147 are plus strand reads.
#'
#' @name import-bam
#' @rdname import-bam
#' @aliases import-bam impBamGAL impBamGAL,character-method
#' @docType methods
#' @param UserDir The relative or full path of folder containing bam files.
#' @param searchArea A \code{GenomicRanges object} that contains the regions of
#' interest
#' @param files use character vector to specify one or more files to import. The
#' default imports all bam files from the directory.
#' @param XStag Setting \code{XStag=TRUE} stores the strand specific
#' information in the mcols slot 'XS'
#' @param verbose makes the function more talkative.
#' @param ... arguments to pass on
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
#' #all files in directory
#' reads <- impBamGAL(pathToFiles,searchArea,verbose=FALSE)
#' #specified files in directory
#' reads <- impBamGAL(pathToFiles,searchArea,
#'				files=c("ERR009160.bam", "ERR009167.bam"),verbose=FALSE)
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
	function(UserDir, searchArea, files = NULL, XStag = FALSE, verbose = TRUE, ...) {

	UserDir <- sub("/$", "", UserDir)

    # Set parameters
    which <- searchArea  #A GRanges, IntegerRangesList, or missing object, from which a IRangesList instance will be constructed.
    what <- scanBamWhat()  #A character vector naming the fields to return. scanBamWhat() returns a vector of available fields. Fields are described on the scanBam help page.
    flag <- scanBamFlag(isUnmappedQuery = FALSE)

    if (XStag) {
        param <- ScanBamParam(flag = flag, which = which, what = what, tag = "XS")  #store ScanBamParam in param.
    } else {
        param <- ScanBamParam(flag = flag, which = which, what = what)  #store ScanBamParam in param.\t\t
    }
    # Point to correct directory and create a BamFileList object
    bamDir <- normalizePath(UserDir)  #Point to the directory containing your Bam files and its respective bam.bai files.

	if(is.null(files)){
		allFiles <- list.files(bamDir, full.names = TRUE)  #list files in a folder.
		bamFiles <- allFiles[grep(".bam$", allFiles)]  #list only the files ending in .bam .
	}else{
		bamFiles <- .mergeDirAndFilename(bamDir, files)
	}

    if (length(bamFiles) == 0)
        stop(paste("No bam files found in", bamDir))
    if (!all(file.exists(paste(bamFiles, ".bai", sep = "")))) {
        if (verbose)
            cat(paste("The bam files in UserDir are required to also have .bam.bai index files in the same directory. Trying to run indexBam function on each"),
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


	UserDir <- sub("/$", "", UserDir)

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
#' @aliases getAlleleCounts getAlleleCounts,GAlignmentsList-method
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


    if (!(is(BamList, "GAlignments") || is(BamList, "GAlignmentsList"))) {
        stop("BamList has to be a GAlignments or GAlignmnetsList object\n")
    }
    # if just one element of, make list (which is a convenient way of
	# handling this input type)
    #
    if (is(BamList, "GAlignments")) {
        BamList <- GAlignmentsList(BamList)
    }

    # check for strand name
    if (!is.character(strand)) {
        stop("strand has to be a character vector")
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
    if (is(GRvariants, "GRangesList")) {
        GRvariants <- unique(unlist(GRvariants, use.names = FALSE))
    }

	#if BamList is not list, make it a list
	if(is(BamList, "GAlignments")){
		BamList <- GAlignmentsList(BamList)
	}

	#Drop seqlevels in BamList that are not in GRvariants
	#seqlevels(BamList,pruning.mode="coarse") <- seqlevels(GRvariants)
	seqinfo(GRvariants) <- merge(seqinfo(GRvariants), seqinfo(BamList))
	seqlevels(GRvariants) <- seqlevelsInUse(GRvariants)


    # check that seqlevels are the same
   # if (!identical(seqlevels(BamList), seqlevels(GRvariants))) {
   #     stop("!identical(seqlevels(BamList), seqlevels(GRvariants))\n")
   # }

    # checking that GRvariants is ok
    if (!is(GRvariants, "GRanges"))
        stop(paste("GRvariants must be a GRanges object, not a",
				   class(GRvariants)))
    if (length(GRvariants) == 0)
        stop("GRvariants was given as an empty GRanges object.",
			 " There can be no Snps retrieved by getAlleleCount then")
    if (any(width(GRvariants) != 1))
        stop("GRvariants can contain only entries of width=1,",
			 " corresponding to SNPs.")

    # checking that verbose is ok
    if (!is.logical(verbose))
        stop(paste("verbose must be a logical, not a", class(verbose)))
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
			if(is.integer(mat)){
				mat <- t(as.matrix(mat))
			}
            if (is.numeric(mat)) {
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

	if(is.list(alleleCountList)){
		MapBiasExpMean <- getDefaultMapBiasExpMean(alleleCountList)
		# make 3D array
		MapBiasExpMean3D <- array(NA, c(length(alleleCountList),
										length(unlist(unique(lapply(alleleCountList,
			rownames)))), 4))  #empty array
		for (i in 1:length(unlist(unique(lapply(alleleCountList, rownames))))) {
			MapBiasExpMean3D[, i, ] <- MapBiasExpMean
		}
	}else if(is.array(alleleCountList)){
		# make 3D array
		MapBiasExpMean3D <- alleleCountList
		mapbiasmat <- t(apply(apply(alleleCountList,c(1,3),sum),
					1, function(x){rank(x,ties.method="first")}))
		mapbiasmat[mapbiasmat==1] <- 0.5
		mapbiasmat[mapbiasmat==2] <- 0.5
		mapbiasmat[mapbiasmat==3] <- 0
		mapbiasmat[mapbiasmat==4] <- 0

		MapBiasExpMean3D <- aperm(array(mapbiasmat, dim=dim(alleleCountList)[c(1,3,2)]),c(1,3,2))

	}
    MapBiasExpMean3D
})

#' Get Gene Area
#'
#' Given a character vector with genesymbols and an OrgDb object, this function
#' returns a GRanges giving the coordinates of the genes.
#'
#' This function is a convenience function that can be used to determine which
#' genomic coordinates to specify to e.g. \code{impBamGAL} when retrieving
#' reads.
#'
#' The function cannot handle genes that do not exist in the annotation. To
#' remove these please set the na.rm=TRUE.
#'
#' @name getAreaFromGeneNames
#' @rdname getAreaFromGeneNames
#' @aliases getAreaFromGeneNames
#' getAreaFromGeneNames,character-method
#' @docType methods
#' @param genesymbols A character vector that contains genesymbols of genes
#' from which we wish to retrieve the coordinates
#' @param OrgDb An \code{OrgDb} object containing gene annotation
#' @param leftFlank A \code{integer} specifying number of additional
#' nucleotides before the genes
#' @param rightFlank A \code{integer} specifying number of additional
#' nucleotides after the genes
#' @param na.rm A \code{boolean} removing genes that returned NA from the
#' annotation
#' @param verbose Setting \code{verbose=TRUE} makes function more talkative
#' @param ... arguments to pass on
#' @return \code{getAreaFromGeneNames} returns a GRanges object with genomic
#' coordinates around the specified genes
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords genes locations
#' @examples
#'
#' #load example data
#' data(ASEset)
#'
#' #get counts at the three positions specified in GRvariants
#' library(org.Hs.eg.db )
#' searchArea<-getAreaFromGeneNames(c('PAX8','TLR7'), org.Hs.eg.db)
#'
NULL

#' @rdname getAreaFromGeneNames
#' @export
setGeneric("getAreaFromGeneNames", function(genesymbols, ...
	){
    standardGeneric("getAreaFromGeneNames")
})

#' @rdname getAreaFromGeneNames
#' @export
setMethod("getAreaFromGeneNames", signature(genesymbols = "character"),
function(genesymbols, OrgDb, leftFlank = 0, rightFlank = 0,
    na.rm = FALSE, verbose = TRUE) {

    # start up sets
    if (!is.character(genesymbols))
        stop(paste("genesymbols must be a character vector, not a", class(genesymbols)))

    if (!is(OrgDb, "OrgDb"))
        stop(paste("OrgDb must be an OrgDb object, not a", class(OrgDb)))

    if (!is.numeric(leftFlank))
        stop(paste("leftFlank must be a numeric vector, not a", class(leftFlank)))
    if (length(leftFlank) != 1)
        stop(paste("leftFlank must be of length 1, not", length(leftFlank)))
    if (leftFlank < 0)
        stop(paste("leftFlank must be equal to or larger than 0"))

    if (!is.numeric(rightFlank))
        stop(paste("rightFlank must be of a numeric vector, not a", class(rightFlank)))
    if (length(rightFlank) != 1)
        stop(paste("rightFlank must be of length 1, not", length(rightFlank)))
    if (rightFlank < 0)
        stop(paste("rightFlank must be equal to or larger than 0"))

    if (!is.logical(verbose))
        stop(paste("verbose must be a logical, not a", class(verbose)))
    if (length(verbose) != 1)
        stop(paste("verbose must be of length 1, not", length(verbose)))

    # retrieving data
    colsFilter <- c("CHR", "CHRLOC", "CHRLOCEND", "SYMBOL")
    s <- suppressWarnings(select(OrgDb, keys = genesymbols, columns = colsFilter, keytype = "SYMBOL"))

    missing <- genesymbols[!genesymbols %in% s[, "SYMBOL"]]
    if (verbose & length(missing) > 0) {
        cat(paste("Did not find information on these", length(missing), "genes:",
            paste(missing, collapse = ", ")), "\n")
    } else {
        if (verbose)
            cat("Found all requested genes in annotation", "\n")
    }

    # remove NAs
    if (na.rm == TRUE) {
        s <- s[!(is.na(s[, "CHR"]) | is.na(s[, "CHRLOC"]) | is.na(s[, "CHRLOCEND"])),
            ]
    } else {
        warningGenes <- s[is.na(s[, "CHR"]) | is.na(s[, "CHRLOC"]) | is.na(s[, "CHRLOCEND"]),
            ]
        if (!nrow(warningGenes) == 0) {
            warning(paste(warningGenes[, "SYMBOL"], "had NAs", "\n"), "you better remove these genes from your 'genesymbols'")
        }
    }

    strand <- rep("+", nrow(s))
    TFstrand <- s[, "CHRLOC"] < 0
    strand[TFstrand] <- "-"

    searchArea <- GRanges(seqnames = paste("chr", s[, "CHR"], sep = ""), ranges = IRanges(abs(s[,
        "CHRLOC"]) - leftFlank, abs(s[, "CHRLOCEND"]) + rightFlank), strand = strand,
        symbol = s[["SYMBOL"]])

    searchArea <- reduce(searchArea, with.revmap = TRUE)
    l <- lapply(searchArea$revmap, function(x) {
        paste(unique(s[x, "SYMBOL"]), collapse = ",")
    })
    mcols(searchArea)[, "symbol"] <- unlist(l)
    # remove the mapping column
    mcols(searchArea)[["revmap"]] <- NULL

    return(searchArea)
})

#' Get rsIDs from locations of SNP
#'
#' Given a GRanges object of SNPs and a SNPlocs annotation, this function
#' attempts to replace the names of the GRanges object entries with rs-IDs.
#'
#' This function is used to try to identify the rs-IDs of SNPs in a GRanges
#' object.
#'
#' @name getSnpIdFromLocation
#' @rdname getSnpIdFromLocation
#' @aliases getSnpIdFromLocation
#' getSnpIdFromLocation,GRanges-method
#' @docType methods
#' @param GR A \code{GRanges} that contains positions of SNPs to look up
#' @param SNPloc A \code{SNPlocs object} containing information on SNP
#' locations (e.g. SNPlocs.Hsapiens.dbSNP.xxxxxxxx)
#' @param return.vector Setting \code{return.vector=TRUE} returns vector with
#' rsIds
#' @param verbose Setting \code{verbose=TRUE} makes function more talkative
#' @param ... arguments to pass on
#' @return \code{getSnpIdFromLocation} returns the same GRanges object it was
#' given with, but with updated with rs.id information.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords SNP rs-id
#' @examples
#'
#' is_32bit_windows <- .Platform$OS.type == "windows" &&
#'                   .Platform$r_arch == "i386"
#' if (!is_32bit_windows && require(SNPlocs.Hsapiens.dbSNP144.GRCh37)) {
#' 	#load example data
#' 	data(ASEset)
#'
#'   #get counts at the three positions specified in GRvariants
#'   updatedGRanges <- getSnpIdFromLocation(rowRanges(ASEset),
#'     SNPlocs.Hsapiens.dbSNP144.GRCh37)
#    rowRanges(ASEset) <- updatedGRanges
#' }
#'
#'
NULL

#' @rdname getSnpIdFromLocation
#' @export
setGeneric("getSnpIdFromLocation", function(GR, ...
	){
    standardGeneric("getSnpIdFromLocation")
})

#' @rdname getSnpIdFromLocation
#' @export
setMethod("getSnpIdFromLocation", signature(GR = "GRanges"),
function(GR, SNPloc, return.vector = FALSE, verbose = TRUE) {

	#genome <- genome(GR)

    if (!is(GR, "GRanges"))
        stop(paste("GR must be a GRanges object, not a", class(GR)))
    if (!is(SNPloc, "SNPlocs"))
        stop(paste("SNPlocs must be a SNPlocs object, not a", class(SNPloc)))
    if (!exists("getSNPlocs"))
        stop("There must exist a function called getSNPlocs, available from the SNPlocs.Hsapiens.dbSNP.xxxxxxx package. Try to load such package.")

	if (any(!(seqlevels(GR) %in% seqlevels(SNPloc)) )){
		stop("one or more seqlevels in GR are not present in SNPloc")
	}
    # remove chr from seqnames if present
   # if (any(grepl("^chr", seqnames(GR))) != length(GR)) {
   #     if (length(grep("^chr", seqnames(GR))) != 0)
   #         stop("seqnames must all begin with 'chr'. In the GR it seemed that some, but not all, seqnames began with chr. Please correct manually")
   #     seqlevels(GR) <- paste("chr", seqlevels(GR), sep = "")
   #     # seqnames(GR)<-seqnames
   # }

    # changing chr to ch to adapt to SNPloc (not anymore)
    #if (length(grep("^chr", seqnames(GR))) == length(GR)) {
    #    # seqnames<-seqnames(GR)
    #    seqlevels(GR) <- sub("^chr", "ch", seqlevels(GR))
    #    # seqlevels(GR)<-as.character(unique(seqnames)) seqlevels(GR) <- levels(seqnames)
    #    # seqnames(GR)<-seqnames
    #}


    #SNPlocThisChr <- getSNPlocs(seqlevels(GR), as.GRanges = TRUE, caching = FALSE)
    #
    #seqlevels(SNPlocThisChr, pruning.mode="coarse") <- seqlevels(GR)
    ## seqlengths(GR) <- seqlengths(SNPlocThisChr)
	#genome(SNPlocThisChr) <- genome(GR)
    #
    #overlaps <- findOverlaps(GR, SNPlocThisChr)
    #
    #if (verbose)
    #    cat(paste("Replacing position-based SNP name with rs-ID for", length(overlaps),
    #        "SNP(s)"), "\n")
    #
    ## replace name in GR for(i in 1:length(overlaps)){
    ## snp<-paste('rs',mcols(SNPlocThisChr[subjectHits(overlaps[i])])[,'RefSNP_id'],sep='')
    ## names(GR)[queryHits(overlaps[i])] <-snp }
    #
    #snp <- paste("rs", mcols(SNPlocThisChr[subjectHits(overlaps)])[, "RefSNP_id"],
    #    sep = "")

	if(any(is.na(genome(GR)))| any(is.null(genome(GR)))){
		genome(GR) <- "NoSetGenome"
	}

	if(!(unique(genome(GR))==unique(genome(SNPloc)))){
		message("setting the GR genome to the same as SNPloc,
				which is a requirement to make overlap calculations")
		genome(GR) <- genome(SNPloc)
	}
	overlap1 <- snpsByOverlaps(SNPloc, GR, type="equal")
	overlap1 <- keepSeqlevels(overlap1, seqlevels(GR))
#	seqinfo(overlap1) <- seqinfo(GR)
#	strand(overlap1) <- "*"
#	mcols(overlap1) <- NULL
#	names(GR) <- NULL

	overlap2 <- findOverlaps(overlap1, GR, type="equal")
    names(GR)[subjectHits(overlap2)] <- mcols(overlap1)[["RefSNP_id"]][queryHits(overlap2)]

    # change back to chr from ch (not needed anymore)
    #seqlevels(GR) <- sub("^ch", "chr", seqlevels(GR))


    if (return.vector) {
        names(GR)
    } else {
        return(GR)
    }

})

#' coverage matrix of GAlignmentsList
#'
#' Get coverage per nucleotide for reads covering a region
#'
#' a convenience function to get the coverage from a list of reads stored in
#' GAlignmnetsList, and returns by default a list with one matrix, and
#' information about the genomic start and stop positions.
#'
#' @name coverageMatrixListFromGAL
#' @rdname coverageMatrixListFromGAL
#' @aliases coverageMatrixListFromGAL
#' coverageMatrixListFromGAL,GAlignmentsList-method
#' @docType methods
#' @param BamList GAlignmentsList containing reads over the region to calculate
#' coverage
#' @param strand strand has to be '+' or '-'
#' @param ignore.empty.bam.row argument not in use atm
#' @param ... arguments to pass on
#' @author Jesper R. Gadin
#' @keywords coverage
#' @examples
#'
#' r <- reads
#' seqlevels(r) <- '17'
#' covMatList <- coverageMatrixListFromGAL(BamList=r, strand='+')
#'
NULL

#' @rdname coverageMatrixListFromGAL
#' @export
setGeneric("coverageMatrixListFromGAL", function(BamList, ...
	){
    standardGeneric("coverageMatrixListFromGAL")
})

#' @rdname coverageMatrixListFromGAL
#' @export
setMethod("coverageMatrixListFromGAL", signature(BamList = "GAlignmentsList"),
function(BamList, strand = "*", ignore.empty.bam.row = TRUE) {

    # If having common start and end points for all gviz track objects the matrix
    # will start on the specific start regardless if there are reads in the bamList
    # or not.

    # TODO, to conveniently access data without loading into memory, the bam file
    # should be read again by using the argument BamPath.

    GAL <- BamList

    # Could be good with a check that the matrix is not longer than CNTNAP2 -
    # 2300000bp long which is the longest gene But will wait with that.
    pstrand = FALSE
    mstrand = FALSE

    if (!is.null(strand)) {
        if (strand == "+") {
            pstrand = TRUE
        } else if (strand == "-") {
            mstrand = TRUE
        } else if (strand == "*") {
            mstrand = TRUE
            pstrand = TRUE
        } else if (strand == "both") {
            mstrand = TRUE
            pstrand = TRUE
        } else {
            stop("strand has to be '+' or '-' if not NULL\n")
        }
    }

    if (!length(seqlevels(GAL)) == 1) {
        stop("can only be one seq level\n")
    }

    # get start and end before filtering on strand, will make things easier
    # downstream.
    suppressWarnings(bamStart <- min(min(start(GAL))))
    suppressWarnings(bamEnd <- max(max(end(GAL))))
    bamWidth <- bamEnd - bamStart + 1

    # if(is.null(start) | is.null(end)){
    start <- bamStart
    end <- bamEnd
    width <- bamWidth
    # }

    if (pstrand) {
        GALp <- GAL[strand(GAL) == "+"]
    }
    if (mstrand) {
        GALm <- GAL[strand(GAL) == "-"]
    }

    if (pstrand) {
        matP <- matrix(0, ncol = (width), nrow = length(GAL))
    }
    if (mstrand) {
        matM <- matrix(0, ncol = (width), nrow = length(GAL))
    }
    if (pstrand) {
        rownames(matP) <- names(GAL)
    }
    if (mstrand) {
        rownames(matM) <- names(GAL)
    }
    ##################################################

    covVecFromGA <- function(GA) {
        mcols(GA) <- NULL
        one <- unlist(grglist(GA))
        covRle <- coverage(one)[[1]]
        cov <- as.integer(window(covRle, start, end))
        cov
    }

    if (pstrand) {
        for (i in 1:length(GALp)) {
            matP[i, ] <- covVecFromGA(GALp[[i]])
        }
    }
    if (mstrand) {
        for (i in 1:length(GALm)) {
            matM[i, ] <- covVecFromGA(GALm[[i]])
        }
    }

    # make mat from matP or matM
    if (strand=="+") {
        mat <- matP
    }
    if (strand=="-") {
        mat <- matM
    }
    if (strand=="*") {
        mat <- matM+matP
    }

    if (strand=="both") {
        mat <- list(plus=matP,minus=matM)
    }

    # store in a list
    if (!is.null(strand)) {
        retList <- list(mat, start, end)
    } else {
        stop("strand must be present")
    }
    # set name on list
    if (!is.null(strand)) {
        names(retList) <- c("mat", "start", "end")
    } else {
        stop("strand must be present")
    }

    retList
})

#' alleleCounts from bam file
#'
#' count alleles before creating ASEse.
#'
#' counts the alleles in a bam file based on GRanges positions.
#'
#' Important excerpt from the details section of the internal applyPileups
#'  function: Regardless of 'param' values, the algorithm follows samtools by
#'  excluding reads flagged as unmapped, secondary, duplicate, or
#'  failing quality control.
#'
#' @name countAllelesFromBam
#' @rdname countAllelesFromBam
#' @aliases countAllelesFromBam
#' countAllelesFromBam,GRanges-method
#' @docType methods
#' @param gr GRanges that contains SNPs of interest
#' @param pathToDir path to directory of bam files
#' @param flag specify one flag to use as filter, default is no filtering.
#' allowed flags are 99, 147, 83 and 163
#' @param scanBamFlag set a custom flag to use as filter
#' @param return.class type of class for the returned object
#' @param verbose makes funciton more talkative
#' @param ... arguments to pass on
#' @author Jesper R. Gadin
#' @keywords allelecount counting
#' @examples
#'
#' data(GRvariants)
#' gr <- GRvariants
#'
#' ##not run at the moment
#' #pathToDir <- system.file('inst/extdata/ERP000101_subset', package='AllelicImbalance')
#' #ar <- countAllelesFromBam(gr, pathToDir)
#'
NULL

#' @rdname countAllelesFromBam
#' @export
setGeneric("countAllelesFromBam", function(gr, ...
	){
    standardGeneric("countAllelesFromBam")
})

#' @rdname countAllelesFromBam
#' @export
setMethod("countAllelesFromBam", signature(gr = "GRanges"),
function(gr, pathToDir, flag=NULL, scanBamFlag=NULL, return.class="array", verbose=TRUE, ...) {

	bamDir <- normalizePath(pathToDir)
	allFiles <- list.files(bamDir, full.names = TRUE)
	bamFiles <- allFiles[grep(".bam$", allFiles)]
	if (length(bamFiles) == 0) {
		stop(paste("No bam files found in", bamDir))
	}
	if (!all(file.exists(paste(bamFiles, ".bai", sep = "")))) {
		if (verbose) {
			cat(paste("The bam files in UserDir are required to also have", ".bam.bai index files.",
				" Trying to run indexBam function on each", "\n"), )
		}
		indexBam(bamFiles)
		if (!all(file.exists(paste(bamFiles, ".bai", sep = "")))) {
			stop("The bam files in UserDir are required to also have", ".bam.bai index files.")
		} else {
			if (verbose) {
				cat(paste("Succesfully indexed all bamFiles in UserDir", pathToDir,
				  "\n"))
			}
		}
	}

	#scanBamFlag
	if(is.null(scanBamFlag) & is.null(flag)){
		flag <- scanBamFlag()
	}
	if(!is.null(scanBamFlag) ){
		if(length(scanBamFlag)==2){
			flag <- scanBamFlag
		}else{
			stop("scanBamFlag has to be the return values from scanBamFlag()")
		}
	}
	#flags can be 99 147 83 or 163
	if(length(flag)==1 & is.null(scanBamFlag)){
		if(! (flag%in%c(99,147,83,163))){
			stop(paste("flag values can only be 99 147 83 or 163",
					   "the input flag value was", flag,sep=" "))
		}

		if(flag==99){

			flag=scanBamFlag( isPaired = TRUE,
					isProperPair = TRUE,
					isFirstMateRead = TRUE,
					isMateMinusStrand = TRUE
					)
		}else if(flag==83){

			flag=scanBamFlag(
					isPaired = TRUE,
					isProperPair = TRUE,
					isFirstMateRead = TRUE,
					isMinusStrand = TRUE
					)
		}else if(flag==147){

			flag=scanBamFlag(
					isPaired = TRUE,
					isProperPair = TRUE,
					isFirstMateRead = FALSE,
					isMinusStrand = TRUE
					)
		}else if(flag==163){

			flag=scanBamFlag(
					isPaired = TRUE,
					isProperPair = TRUE,
					isFirstMateRead = FALSE,
					isMateMinusStrand = TRUE
					)
		}
	}


	if(verbose){
		cat("sam flag used:",flag,"\n")
	}
	fls <- PileupFiles(bamFiles)

	countF <-
		function(x){
		x[["seq"]][-5,,1]

	}

	which <- gr
	p1 <- ApplyPileupsParam(flag=flag,
							which=which,
							minBaseQuality = 0L,
							what="seq",
							yieldBy = "position",
							yieldAll=TRUE,
							maxDepth=.Machine$integer.max,
							...
							)

	res <- applyPileups(fls, countF, param=p1)

	ar <- array(unlist(res), dim=c(4, length(fls), length(which)),
				dimnames=list(c("A","C","G","T"),
							  names(fls),
							  names(which)))
	ar <- aperm(ar,dim=c(3,2,1))

	ar
})

#' ASEset from bam file
#'
#' count alleles and create an ASEset direct from bam file instead of reading into R first.
#'
#' counts the alleles in a bam file based on GRanges positions.
#'
#'
#' @name ASEsetFromBam
#' @rdname ASEsetFromBam
#' @aliases ASEsetFromBam
#' ASEsetFromBam,GRanges-method
#' @docType methods
#' @param gr GenomicRanges of SNPs to create ASEset for
#' @param PE if paired end or not (default: TRUE)
#' @param pathToDir Directory of bam files with index in same directory
#' @param strandUnknown default: FALSE
#' @param ... passed on to ASEsetFromBam function
#' @param flagsMinusStrand flags that mark reads coming from minus strand
#' @param flagsPlusStrand flags that mark reads coming from plus strand
#' @author Jesper R. Gadin
#' @keywords ASEset
#' @examples
#'
#' data(GRvariants)
#' gr <- GRvariants
#'
#' ##no execution at the moment
#' #pathToDir <- system.file('inst/extdata/ERP000101_subset', package='AllelicImbalance')
#' #a <- ASEsetFromBam(gr, pathToDir)
#'
NULL

#' @rdname ASEsetFromBam
#' @export
setGeneric("ASEsetFromBam", function(gr, ...
	){
    standardGeneric("ASEsetFromBam")
})

#' @rdname ASEsetFromBam
#' @export
setMethod("ASEsetFromBam", signature(gr = "GRanges"),
 function(gr, pathToDir,PE=TRUE, flagsMinusStrand=c(83,163), flagsPlusStrand=c(99,147), strandUnknown=FALSE, ...) {

	if(!PE){
		stop("no support for SE atm")
	}

	if(PE==TRUE){
		#minus strand
		arm1 <- countAllelesFromBam(gr, pathToDir, flag=83)
		arm2 <- countAllelesFromBam(gr, pathToDir, flag=163)
		arm <- arm1 + arm2

		#plus strand
		arp1 <- countAllelesFromBam(gr, pathToDir, flag=99)
		arp2 <- countAllelesFromBam(gr, pathToDir, flag=147)
		arp <- arp1 + arp2
	}

	#ASEsetFromArray
	if(!strandUnknown){
		a <- ASEsetFromArrays(gr, countsPlus = arp,
			countsMinus = arm)

	}else{
		a <- ASEsetFromArrays(gr, countsUnknown = arp+arm)
	}
	a
})

#' makes masked fasta reference
#'
#' Replaces all selected positions in a fasta file with the character N
#'
#' @name makeMaskedFasta
#' @rdname makeMaskedFasta
#' @aliases makeMaskedFasta
#' makeMaskedFasta,character-method
#' @docType methods
#' @param fastaIn character string of the path for the fasta file to be used
#' @param fastaOut character string of the path for the masked fasta file (no extension)
#' @param posToReplace GRanges object with the genomic ranges to replace
#' @param splitOnSeqlevels write on file for each seqlevel to save memory
#' @param verbose makes function more talkative
#' @param ... arguments to pass on
#' @author Jesper R. Gadin
#' @keywords masked fasta reference
#' @examples
#'
#' data(ASEset.sim)
#' gr <- rowRanges(ASEset.sim)
#' fastaIn <- system.file('extdata/hg19.chr17.subset.fa', package='AllelicImbalance')
#' makeMaskedFasta(fastaIn=fastaIn, fastaOut="fastaOut",posToReplace=gr)
#'
#'
NULL

#' @rdname makeMaskedFasta
#' @export
setGeneric("makeMaskedFasta", function(fastaIn, ...
	){
    standardGeneric("makeMaskedFasta")
})

#' @rdname makeMaskedFasta
#' @export
setMethod("makeMaskedFasta", signature(fastaIn = "character"),
function(fastaIn, fastaOut, posToReplace, splitOnSeqlevels=TRUE, verbose=TRUE){

	#does the inFasta file exist?
	if(!file.exists(fastaIn)){
		stop("fasta infile doesnt exist")
	}

	#make FaFile
	fl <- FaFile(fastaIn)

	#check if the index file is present otherwise tell the user to use the
	#indexFa(FaFile("pathToReference")) command
	if(!file.exists(paste(fastaIn,".fai",sep=""))){
		cat("could not find index file\n")
		cat("creates a new index file")
		indexFa(FaFile(fastaIn))
		cat("finished creating new index file")
	}
	#indexFa(fl) #creates a new file as index in the same directory but
	#with extension *.fai

	#open,scan,close file
	open(fl)
	#check if seqleveles are present
	fa.info <- scanFaIndex(fl)
	if(!(all(seqlevels(posToReplace) %in% seqlevels(fa.info)))){
		close(fl)
		stop("seqlevels in object x are not in fasta index file")
	}

	###############LOOOP
	chrs <- seqlevels(posToReplace)
	for (chr in chrs){
	#	chr <- "chr1"

		searchArea <- fa.info[seqnames(fa.info)==chr]
		seq <- scanFa(fl,param=searchArea)

		#replace the SNPs with N
		toReplace <- unique(as.integer(ranges(posToReplace)))
		seq <- seq[[1]]
		seq[toReplace] <- "N"

		if(verbose){cat("replaced", length(toReplace),"instances with N\n")}

		#write new file
		#library("seqinr")
		outfile <- paste(fastaOut,chr,".fa",sep="")

		write.fasta(as.character(seq),names=chr, file.out=outfile, nbchar=80)
		if(verbose){cat("wrote chr",chr,"to file\n")}

		gc()
	}
	close(fl)
	if(verbose){cat("all chromosomes written to file\n")}
})


#' global analysis wrapper
#'
#' A wrapper to make a global analysis based on paths for BAM, VCF and GFF files
#'
#' @name gba
#' @rdname gba
#' @aliases gba
#' gba,character-method
#' @docType methods
#' @param pathBam path to bam file
#' @param pathVcf path to vcf file
#' @param pathGFF path to gff file
#' @param verbose makes function more talkative
#' @param ... arguments to pass on
#' @author Jesper R. Gadin
#' @keywords global wrapper
#' @examples
#'
#' #empty as function doesn't exist
#'
NULL

#' @rdname gba
#' @export
setGeneric("gba", function(pathBam, ...
	){
    standardGeneric("gba")
})

#' @rdname gba
#' @export
setMethod("gba", signature(pathBam = "character"),
function(pathBam,pathVcf,pathGFF=NULL, verbose){

	#summarize counts

	#detectAI

})


#' snp quality data
#'
#' Given the positions of known SNPs, this function returns allele quality from
#' a BamGRL object
#'
#' This function is used to retrieve the allele quality strings from specified positions
#' in a set of RNA-seq reads. The \code{BamList} argument will typically have
#' been created using the \code{impBamGAL} function on bam-files. The
#' \code{GRvariants} is either a GRanges with user-specified locations or else
#' it is generated through scanning the same bam-files as in \code{BamList} for
#' heterozygote locations (e.g. using \code{scanForHeterozygotes}). The
#' GRvariants will currently only accept locations having width=1,
#' corresponding to bi-allelic SNPs. The strand type information will be kept in the
#' returned object. If the strand is marked as unknown "*", it will be forced to the "+"
#' strand.
#'
#' quaity information is extracted from the BamList object, and requires the presence of
#' mcols(BamList)[["qual"]] to contain quality sequences.
#'
#' @name getAlleleQuality
#' @rdname getAlleleQuality
#' @aliases getAlleleQuality getAlleleQuality,GAlignmentsList-method
#' @docType methods
#' @param BamList A \code{GAlignmentsList object} or \code{GRangesList object}
#' containing data imported from a bam file
#' @param GRvariants A \code{GRanges object} that contains positions of SNPs to
#' retrieve.
#' @param fastq.format default 'illumina.1.8'
#' @param return.class 'list' or 'array'
#' @param verbose Setting \code{verbose=TRUE} makes function more talkative
#' @param ... parameters to pass on
#' @return \code{getAlleleQuality} returns a list of several data.frame objects,
#' each storing the count data for one SNP.
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords allele quality
#' @examples
#'
#' #load example data
#' data(reads)
#' data(GRvariants)
#'
#' #get counts at the three positions specified in GRvariants
#' alleleQualityArray <- getAlleleQuality(BamList=reads,GRvariants)
#'
#' #place in ASEset object
#' alleleCountsArray <- getAlleleCounts(BamList=reads,GRvariants,
#'                      strand='*', return.class="array")
#'
#' 	a <- ASEsetFromArrays(GRvariants, countsUnknown = alleleCountsArray)
#' 	aquals(a) <- alleleQualityArray
NULL

#' @rdname getAlleleQuality
#' @export
setGeneric("getAlleleQuality", function(BamList, ...
	){
    standardGeneric("getAlleleQuality")
})

#' @rdname getAlleleQuality
#' @export
setMethod("getAlleleQuality", signature(BamList = "GAlignmentsList"),
function(BamList, GRvariants, fastq.format = "illumina.1.8",
						return.class = "array", verbose = TRUE, ...) {

	if(!return.class %in% c("array")){
		stop("return.class has to be array")
	}

    if (!(is(BamList, "GAlignments") || is(BamList, "GAlignmentsList"))) {
        stop("BamList has to be a GAlignments or GAlignmnetsList object\n")
    }
    # if just one element of, make list (which is a convenient way of
	# handling this input type)
    #
    if (is(BamList, "GAlignments")) {
        BamList <- GAlignmentsList(BamList)
    }

    # if the user sent in the GRangesList for GRvariants,
	# take out only the unique entries.
    #
    if (is(GRvariants, "GRangesList")) {
        GRvariants <- unique(unlist(GRvariants, use.names = FALSE))
    }

	#if BamList is not list, make it a list
	if(is(BamList, "GAlignments")){
		BamList <- GAlignmentsList(BamList)
	}

	#Drop seqlevels in BamList that are not in GRvariants
	#seqlevels(BamList,pruning.mode="coarse") <- seqlevels(GRvariants)
	seqinfo(GRvariants) <- merge(seqinfo(GRvariants), seqinfo(BamList))
	seqlevels(GRvariants) <- seqlevelsInUse(GRvariants)


    # check that seqlevels are the same
   # if (!identical(seqlevels(BamList), seqlevels(GRvariants))) {
   #     stop("!identical(seqlevels(BamList), seqlevels(GRvariants))\n")
   # }

    # checking that GRvariants is ok
    if (!is(GRvariants, "GRanges"))
        stop(paste("GRvariants must be a GRanges object, not a",
				   class(GRvariants)))
    if (length(GRvariants) == 0)
        stop("GRvariants was given as an empty GRanges object.",
			 " There can be no Snps retrieved by getAlleleCount then")
    if (any(width(GRvariants) != 1))
        stop("GRvariants can contain only entries of width=1,",
			 " corresponding to SNPs.")

    # checking that verbose is ok
    if (!is.logical(verbose))
        stop(paste("verbose must be a logical, not a", class(verbose)))
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

	#choose format
	if(fastq.format=="illumina.1.8"){
		dim3 <- c("!","\"","#","$","%","&","'","(",")","*","+",",","-",".","\\","/","0","1","2","3","4","5","6","7","8","9",":",";","<","=",">","?","@","A","B","C","D","E","F","G","H","I","J")
	}


	#empty array that handles only four nucleotides + one del columns
    dimnames = list(snpNames, names(BamList), dim3, c("+","-"))
    ar1 <- array(NA, c(length(GRvariants), length(BamList), length(dim3), 2),
				 dimnames = dimnames)

    # force all "*" on "+" strand
	un1 <- strand(BamList)=="*"
	un2 <- any(un1)
	if(any(un2)){
		strand(BamList)[un2][un1[un2]] <- "+"
	}

    for (j in 1:length(names(BamList))) {
        sample <- names(BamList)[j]
        if (verbose)
            cat("sample ", sample, "\n")

		my_IGPOI <- GRvariants

		# + strand
        gal <- BamList[[j]][strand(BamList[[j]]) == "+"]
		seqlevels(gal) <- seqlevels(my_IGPOI)

		nuclpiles <- pileLettersAt(mcols(gal)$qual, seqnames(gal), start(gal), cigar(gal),
							                                 my_IGPOI)
        # fill array
        nstr <- factor(unlist(strsplit(as.character(nuclpiles), "")),
					   levels=dim3, labels=dim3)
		levels(nstr) <- dim3

        for (k in 1:length(nuclpiles)) {
			tbl <- table(factor(unlist(strsplit(as.character(nuclpiles[k]), "")),
					   levels=dim3, labels=dim3))
			ar1[k, j, names(tbl), "+"] <- as.integer(tbl)
        }

		# - strand
        gal <- BamList[[j]][strand(BamList[[j]]) == "-"]
		seqlevels(gal) <- seqlevels(my_IGPOI)

		nuclpiles <- pileLettersAt(mcols(gal)$qual, seqnames(gal), start(gal), cigar(gal),
							                                 my_IGPOI)

        # fill array
        nstr <- factor(unlist(strsplit(as.character(nuclpiles), "")),
					   levels=dim3, labels=dim3)
		levels(nstr) <- dim3

        for (k in 1:length(nuclpiles)) {
			tbl <- table(factor(unlist(strsplit(as.character(nuclpiles[k]), "")),
					   levels=dim3, labels=dim3))
			ar1[k, j, names(tbl), "-"] <- as.integer(tbl)
        }
    }

	if (return.class == "array") {
        ar1
    } else {
        cat("return.class unknown\n Nothing will be returned from function!")
    }
})


#' genotype2phase
#'
#' used to convert the genomatrix from the visually friendly matrix to phase array.
#'
#' To not introduce redundant information in the ASEset object, the genotype matrix is
#' translated to a phase matrix, containing the same information.
#' Does not allow tri-allelic or multi-allelic SNPs, and if present the multi-allelic
#' SNPs will lose the least occuring genotype.
#'
#' This function can handle indels, but if the reference allele is not provided, the
#' rank matrix which is temporary created might use lots of memory, depending on the
#' amount of indels among the genotypes. As conclusion, it is preferable to send in
#' reference genome when converting to phase.
#'
#' levels information is only important if the reference allele has to be guessed,
#' and so if reference information is provided, the levels argument can be ignored.
#'
#' @name genotype2phase
#' @rdname genotype2phase
#' @aliases genotype2phase,matrix-method
#' @docType methods
#' @param x matrix see examples
#' @param ref reference alleles
#' @param return.class 'array' or 'list'
#' @param levels vector of expected alleles
#' @param ... pass on additional param
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords phase
#' @examples
#'
#' #load example data
#' data(genomatrix)
#' data(ASEset)
#' p <- genotype2phase(genomatrix, ref(ASEset))
#'
NULL

#' @rdname genotype2phase
#' @export
setGeneric("genotype2phase", function(x, ...
	){
    standardGeneric("genotype2phase")
})

#' @rdname genotype2phase
#' @export
setMethod("genotype2phase", signature(x = "matrix"),
		function(x, ref=NULL, return.class="array", levels=c("A","C","G","T") , ...
	){

	#check x
	if(!is.matrix(x)){
		stop("x is not of a matrix")
	}

	#checks return.class
	if(return.class=="list" & !is.null(ref)){
		stop("reference is already known, no need to use return.class='list'")
	}

	#check ref
	if(!is.null(ref)){
		if(!length(ref)==nrow(x)){
		 stop("reference length is not equal to genotype matrix nrow")
		}
	}

	sgm <- .splitGenotypeMatrix(x)
	if(is.null(ref)){
		#if length of levels is to long the rank matrix risks being very huge
		sgc <- .splitGenotypeCount(sgm, levels)
		sgr <- .splitGenotypeRank(sgc)
		ref <- .splitGenotypePickRefAllele(sgr)
		alt <- .splitGenotypePickAltAllele(sgr, ref)
	}

	p <-  array(c((!sgm == ref)*1, rep(0, length(sgm[,,1]))),dim= c(nrow(x), ncol(x), 3))

	if(return.class=="list"){
		list(phase=p,ref=ref,alt=alt)
	}else if(return.class=="array"){
		p
	}
})

### -------------------------------------------------------------------------
### helpers for genotype2phase
###

#split a genotype matrix into an array of two dimensions(one for each allele)
.splitGenotypeMatrix <-
	function(genomatrix, ...)	{
	str <- unlist(strsplit(genomatrix,"/"))
	inx <- which(is.na(str))
	val <- c(str, rep(NA,length(inx)))
	id  <- sort(c( seq_along(str), inx+0.5), index.return=TRUE)$ix

	aperm(array(val[id], dim=c(2, nrow(genomatrix), ncol(genomatrix))),c(2,3,1))
}

#from a genotype array of two dimensions, count occurences of each allele for
#each snp
.splitGenotypeCount <-
	function(x, levels=c("A", "C", "G", "T"), ...)
{
	mat <- matrix(aperm(x, c(3,2,1)), ncol(x)*2, nrow(x))
	t(apply(mat, 2, function(x){table(factor(x, levels=levels))}))
}

#from matrix rank the alleles for each SNP
.splitGenotypeRank <-
	function(x, levels=colnames(x), ...)
{
	mat <- t(apply(x, 1, function(x){sort(x, decreasing=TRUE, index.return=TRUE)$ix}))
	matrix(levels[mat], nrow(x), ncol(x), dimnames=list(NULL,paste("r",1:length(levels),sep="")))
}

#pick a random ref allele from two most occured
.splitGenotypePickRefAllele <- function(x, ...){
	apply(x[,1:2], 1, sample,1)
}

#pick a random alt allele from two most occured (which is not ref allele)
.splitGenotypePickAltAllele <- function(x, ref, ...){
	x[,1:2][!x[,1:2] == ref]
}


#' phase2genotype
#'
#' Convert the phase from the internally stored phase, ref and alt information
#'
#' To not introduce redundant information in the ASEset object, the genotype matrix is
#' accessed from the phase matrix, which together with ref and alt allele information
#' contains the same information(not taken into account three-allelic or more SNPs).
#'
#' The genotype matrix retrieved from an ASEset object can differ from the genotype matrix
#' stored in the object if reference and alternative alleles were not used or has changed
#' since the phase genotype matrix was stored. Basically, it is preferable to provide
#' reference and alternative information when storing the genotype matrix.
#'
#' If possible, it is better to not use a genotype matrix, but instead relying completely
#' on storing a phase matrix(or array) together with reference and alternative allele
#' information.
#'
#' @name phase2genotype
#' @rdname phase2genotype
#' @aliases phase2genotype,array-method
#' @docType methods
#' @param x array see examples
#' @param ref reference allele vector
#' @param alt alternative allele vector
#' @param return.class 'matrix' or 'array'
#' @param ... pass on additional param
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords phase
#' @examples
#'
#' #load example data
#' data(ASEset)
#' data(genomatrix)
#' p <- genotype2phase(genomatrix, ref(ASEset), return.class="array")
#' ref <- ref(ASEset)
#' alt <- inferAltAllele(ASEset)
#'
#' gt <- phase2genotype(p, ref, alt, return.class="matrix")
#'
NULL

#' @rdname phase2genotype
#' @export
setGeneric("phase2genotype", function(x, ...
	){
    standardGeneric("phase2genotype")
})

#' @rdname phase2genotype
#' @export
setMethod("phase2genotype", signature(x = "array"),
		function(x, ref, alt, return.class="matrix", ...
	){

	gta <- .phaseArray2genotypeArray(x, ref, alt)

	if(return.class=="matrix"){
		.genotypeArray2genotypeMatrix(gta)
	}else if(return.class=="array"){
		gta
	}else{ stop("return.class type doesnt exist")}

})

### -------------------------------------------------------------------------
### helpers for phase2genotype
###
#NAs will be kept
.phaseArray2genotypeArray <- function(x, ref, alt){
	mat <- matrix(alt, nrow(x), ncol(x))
	pat <- matrix(alt, nrow(x), ncol(x))
	pha <- matrix("/", ncol=ncol(x), nrow=nrow(x))

	namat <- is.na(x[,,1])
	napat <- is.na(x[,,2])

	mat[x[,,1]==0 & !namat] <- matrix(ref, nrow(x), ncol(x))[x[,,1]==0 & !namat]
	pat[x[,,2]==0 & !napat] <- matrix(ref, nrow(x), ncol(x))[x[,,2]==0 & !napat]

	mat[namat] <- NA
	pat[napat] <- NA

	pha[x[,,3]==1] <- "|"

	array(c(mat, pat, pha), dim=c(dim(x)),
		  dimnames=list(NULL, NULL, c("mat","pat", "phased")))
}

#NAs will be kept and merged
.genotypeArray2genotypeMatrix <- function(x, ...)
{
	napha <- is.na(x[,,3])
	x[,,3][napha] <- "/"

	gv <- paste(x[,,1], x[,,3], x[,,2], sep="")
	namat <- is.na(x[,,1])
	napat <- is.na(x[,,2])

	gv[namat | napat] <- NA
	matrix(gv, nrow(x), ncol(x))
}

##
## DNAStringSet2character
##
## forces e.g. simplelist of variants to an atomic character vector
##
## @name DNAStringSet2character
## @rdname DNAStringSet2character
## @aliases DNAStringSet2character
## DNAStringSet2character,DNAStringSet-method
## @docType methods
## @param x DNAStringSet
## @param verbose makes function more talkative
## @param ... arguments to pass on
## @author Jesper R. Gadin
## @keywords global wrapper
## @examples
##
## #empty as function doesn't exist
##
#NULL
#
## @rdname gba
## @export
#setGeneric("DNAStringSet2character", function(x, ...
#	){
#    standardGeneric("DNAStringSet2character")
#})
#
## @rdname gba
## @export
#setMethod("DNAStringSet2character", signature(x = "character"),
#function(x, verbose){
#
#	if(any(width(x)>1)){
#		stop("One or more variants are not bi-allelic SNPs")
#	}
#	as.character(x)
#
#)
#
#
#isSnp <- function(x) {
#	refSnp <- nchar(ref(x)) == 1L
#	a <- alt(x)
#	altSnp <- elementNROWS(a) == 1L
#	ai <- unlist(a[altSnp]) # all length 1, so unlisting is 1:1 map
#	altSnp[altSnp] <- nchar(ai) == 1L & (ai %in% c("A", "C", "G", "T"))
#	refSnp & altSnp
#}
#
#
