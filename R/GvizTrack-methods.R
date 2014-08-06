setMethod("ASEDAnnotationTrack",
    signature(x = "ASEset"),
    function (x,
		GR=rowData(x),
		type="fraction",
		strand="nonStranded",
		...
	){

		#change to "*"
		#if(strand=="nonStranded"){strand <- "*"} not possile as long as nonStranded is an option

		#check genome
		if(is.null(genome(x)) | is.na(genome(x))){
			stop(paste("genome have be set for object x", "e.g. genome(x) <- \"hg19\" "))	
		}

		#check seqnames has length=0
		if(!(length(seqlevels(x))==1)){stop("This function can only use objects with one seqlevel")}

		if(sum(strand=="+"| strand=="-")==0){
			stop("strand must be plus or minus at the moment")
		}
		if(!nrow(x)==1){
			
			GR <- GRanges(seqnames=seqlevels(x),ranges=IRanges(start=min(start(x)),end=max(end(x))),strand=strand, genome=genome(x))
		
		}
		
		ranges <- rowData(x)

		colnames(x)<- 1:ncol(x)

		details <- function(identifier, ...) {

		type <- get("type",envir=AllelicImbalance.extra)
		arank <- get("arank",envir=AllelicImbalance.extra)
		afraction <- get("afraction",envir=AllelicImbalance.extra)
		acounts <- get("acounts",envir=AllelicImbalance.extra)

			if(type == "fraction"){
				print(barplot.lattice.fraction(identifier,afraction, arank, ... ), 
				newpage = FALSE,
				prefix = "plot")

			}else if(type == "count"){
				print(barplot.lattice.counts(identifier, arank, acounts, ...), 
				newpage = FALSE,
				prefix = "plot")
			}

		}

		#pick out plot data from ASEset
		#strand="+"
		AllelicImbalance.extra <- new.env(parent = emptyenv())

		assign("acounts", alleleCounts(x,strand=strand), envir = AllelicImbalance.extra)
		assign("arank", arank(x,strand=strand), envir = AllelicImbalance.extra)
		assign("afraction", fraction(x, strand=strand), envir = AllelicImbalance.extra)
		assign("type", type , envir = AllelicImbalance.extra)

		#plot the fraction
		deTrack <- AnnotationTrack(range = ranges, genome = genome(x),
			id = rownames(x), name = paste("Gviz locationplot",type ),
			stacking = "squish", fun = details)
		deTrack
	}
    )


setMethod("CoverageDataTrack",
    signature(x = "ASEset"),
    function (x,
		GR=NULL,
		BamList=NULL,
		strand=NULL,
		verbose=TRUE,
		...
	){
		#GR is not in use atm. Missing is a subset of the return matrix based on the GR values.

		if(!is.null(strand)){
			if(strand=="+"){pstrand=TRUE}
			else if(strand=="-"){mstrand=TRUE}
			else{stop("strand has to be '+' or '-' if not NULL\n")}
		}else{stop("strand has to be '+' or '-' if not NULL\n")}

		#check genome
		if(is.null(genome(x)) | is.na(genome(x))){
			stop(paste("genome have be set for object x", "e.g. genome(x) <- \"hg19\" "))	
		}

		#check seqnames has length=0
		if(!(length(seqlevels(x))==1)){stop("This function can only use objects with one seqlevel")}

		if(!is.null(strand)){
			if(strand=="+"){pstrand=TRUE}
			else if(strand=="-"){mstrand=TRUE}
			else{stop("strand has to be '+' or '-' if not NULL\n")}
		}

		GR <- GRanges(seqnames=seqlevels(x),ranges=IRanges(start=min(start(x)),end=max(end(x))),strand=strand)


#	start <- start(GR)
#	end <- end(GR)
#	chr <- seqnames(GR)	



	if(is.null(BamList)){
		stop("must include GappedAlignmentsList as BamList ")
	}else{
		#check that only one chromosome is present
		if(!length(seqlevels(BamList))==1){stop("can only be one seq level\n")}

		covMatList <- coverageMatrixListFromGAL(BamList,strand)
	}

		trackList <- list() 

		#set matP or matM to only mat
		#if(strand=="+"){mat <- covMatList[["matP"]]}
		#if(strand=="-"){mat <- covMatList[["matM"]]}

		mat <- covMatList[["mat"]]
		start <- covMatList[["start"]]
		end <- covMatList[["end"]]

		#prepare Gviz dtracks
		for(j in 1:nrow(mat) ){
			trackList[[length(trackList)+1]] <- 
				DataTrack(
				  data = mat[j,], 
			  	  start=start:end,
			  	  width=1, 
			  	  chromosome = seqlevels(x), 
			  	  genome = genome(x), 
			  	  name = rownames(mat)[j], 
			  	  type="s"
				)
		}

		trackList
	}
    )




