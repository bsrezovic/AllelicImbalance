setMethod("ASEDAnnotationTrack",
    signature(x = "ASEset"),
    function (x,
		type="fraction",
		strand="nonStranded"
	){

		ASEset <- x
		ranges <- rowData(ASEset)
		#strand <- "+"
		if(strand=="nonStranded"){
			strand(ranges) <- "*" 
		}else{
			strand(ranges) <- strand 
		}
		
		
		colnames(ASEset)<- 1:ncol(ASEset)

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

		assign("acounts", alleleCounts(ASEset,strand=strand), envir = AllelicImbalance.extra)
		assign("arank", arank(ASEset,strand=strand), envir = AllelicImbalance.extra)
		assign("afraction", fraction(ASEset, strand=strand), envir = AllelicImbalance.extra)
		assign("type", type , envir = AllelicImbalance.extra)

		#plot the fraction
		deTrack <- AnnotationTrack(range = ranges, genome = "hg19",
			chromosome = 17, id = rownames(ASEset), name = paste("Gviz locationplot",type ),
			stacking = "squish", fun = details)
		deTrack
	}
    )


setMethod("CoverageDataTrack",
    signature(x = "ASEset"),
    function (x,
		BamList=NULL,
		strand=NULL,
		start=NULL,
		end=NULL,
		chr=NULL,
		genome="hg19",
		verbose=TRUE,
		...
	){

	if(!is.null(strand)){
		if(strand=="+"){pstrand=TRUE}
		else if(strand=="-"){mstrand=TRUE}
		else{stop("strand has to be '+' or '-' if not NULL\n")}
	}
	if(is.null(BamList)){
		stop("must support BamList. In future it will be possible to use coverage matrix stored in ASEset")
	}else{
		#check that only one chromosome is present
		if(!length(seqlevels(BamList))==1){stop("can only be one seq level\n")}

		covMatList <- coverageMatrixListFromGAL(BamList,start,end,strand)
	}

		trackList <- list() 

		#set matP or matM to only mat
		#if(strand=="+"){mat <- covMatList[["matP"]]}
		#if(strand=="-"){mat <- covMatList[["matM"]]}
		mat <- covMatList[["mat"]]

		if(is.null(start) | is.null(end)){
			start <- covMatList[["start"]]
			end <- covMatList[["end"]]
		}

		#prepare Gviz dtracks
		for(j in 1:nrow(mat) ){
			trackList[[length(trackList)+1]] <- 
				DataTrack(
				  data = mat[j,], 
			  	  start=start:end,
			  	  width=1, 
			  	  chromosome = seqlevels(x), 
			  	  genome = genome, 
			  	  name = rownames(mat)[j], 
			  	  type="s"
				)
		}

		trackList
	}
    )




