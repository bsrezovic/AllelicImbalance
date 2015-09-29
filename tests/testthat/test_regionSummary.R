context("test units for method regionSummary")

test_that(paste("checking .aiRegionDirectionFromMaternalPhaseMapBiasAndIndex"), {

	#####################
	# Test 1
	#####################
	#Samples on rows
	#SNPs on cols
	fr <- matrix(c(0.3, 0.3, 0.4, 0.5,
				   0.7, 0.6, 0.5, 0.6,
		  	       NA ,  NA,  NA,  NA), ncol=3, nrow=4, byrow=FALSE)
	#same dim as fr
	mb <- array(0.5,dim=dim(fr))
	idx <- c(1,1,2)

	#prepare expected data 
	e1 <- c(1,1,0,1,0,0,0,0)
	e2 <- c(1,1,1,0,0,0,0,0)
	exp <- array(c(e1, e2),dim=c(4,2,2))

	#run tests
	res <- .aiRegionDirectionFromMaternalPhaseMapBiasAndIndex(fr, mb, idx)
	
	#test equality
    expect_that(exp, equals(res))

})

test_that(paste("checking .regionStatisticsFromMatrixAndIndex"), {

	#####################
	# Test 1 - mean
	#####################
	#Samples on rows
	#SNPs on cols
	mat <- matrix(c(0.3, 0.3, 0.4, 0.5,
				   0.7, 0.6, 0.5, 0.6,
		  	       NA ,  NA,  NA,  NA), ncol=3, nrow=4, byrow=FALSE)
	#same dim as fr
	idx <- c(1,1,2)

	#prepare expected data 
	e1 <- c(0.5, 0.45, 0.45, 0.55, NA, NA, NA, NA)
	exp <- matrix(e1, nrow=4, ncol=2)
	dimnames(exp) <- list(1:4,1:2)

	#run tests
	res <- .regionStatisticsFromMatrixAndIndex(mat, idx, mean)
	
	#test equality
    expect_equal(exp, res, tolerance=1e-6)

})

test_that(paste("checking .deltaFromFractionMatrixAndMapBias"), {

	#####################
	# Test 1 - mean
	#####################
	#Samples on rows
	#SNPs on cols
	ar <- matrix(c(0.3, 0.3, 0.4, 0.5,
				   0.7, 0.6, 0.5, 0.6,
		  	       NA ,  NA,  NA,  NA), ncol=3, nrow=4, byrow=FALSE)
	mb <- array(0.5, dim=dim(ar))

	#prepeare expected data
	e1 <- c(0.2, 0.2, 0.1, 0.0, 0.2, 0.1, 0.0, 0.1, NA, NA, NA, NA)
	exp <- matrix(e1, nrow=4, ncol=3)

	#run tests
	res <- .deltaFromFractionMatrixAndMapBias(ar, mb)
	
	#test equality
    expect_that(exp, equals(res))

})

test_that(paste("checking .selectRegionAndTransferIndexToASEset"), {

	#####################
	# Prepare data 1 
	#####################
	data(ASEset) 
	a <- ASEset

	#####################
	# Test 1 - one level region
	#####################
	reg <- c(split(granges(a)[c(1,2,2,3,1,1)], c(1,1,2,2,3,3)))
			 
	start(reg)[[3]] <- c(1,2)
	end(reg)[[3]] <- c(1,2)

	#prepeare expected data
	reg <- .unlistGRangesListAndIndex(reg)

	#prepeare expected data
	exp <- a[c(1,2,2,3)]
	idx <- c(1,1,2,2)
	idn <- paste(c(1,1,2,2))
	mcols(exp)[["ASEsetIndex"]] <- c(1,2,2,3)
	mcols(exp)[["regionIndex"]] <- DataFrame(lvl1=idx)
	mcols(exp)[["regionIndexName"]] <- DataFrame(lvl1=idn)

	#run tests
	res <- .selectRegionAndTransferIndexToASEset(a, reg)
	
	#test equality
    expect_that(exp, equals(res))

})


test_that(paste("checking .unlistGRangesListAndIndex"), {

	#####################
	# Prepare data 1 
	#####################
	data(ASEset) 
	a <- ASEset

	#####################
	# Test 1 
	#####################
	
	#add one region that is missing
	reg <- c(split(granges(a)[c(1,2,2,3,1,1)], c(1,1,2,2,3,3)),
			 split(granges(a)[c(1,2,2,3,1,1)], c(1,1,2,2,3,3)))
	start(reg)[[3]] <- c(1,2)
	end(reg)[[3]] <- c(1,2)

	#prepeare expected data
	idx <- c(1,1,2,2,3,3,4,4,5,5,6,6)
	idn <- paste(c(1,1,2,2,3,3,1,1,2,2,3,3))
	exp <- unlist(reg)
	mcols(exp)[["regionIndex"]] <- DataFrame(lvl1=idx)
	mcols(exp)[["regionIndexName"]] <- DataFrame(lvl1=idn)

	#run tests
	res <- .unlistGRangesListAndIndex(reg)
	
	#test equality
    expect_that(exp, equals(res))

})

test_that(paste("checking .makeRegionGRangesFromASEsetWithRegionIndex"), {

	#####################
	# Prepare data 1 
	#####################
	data(ASEset) 
	a <- ASEset

	#####################
	# Test 1 
	#####################
	
	#add one region that is missing
	x <- a[c(1,2,3,2,2,2,3,2,1,3)]
	mcols(x)[["ASEsetIndex"]] <- c(1,2,3,2,2,2,3,2,1,3)
	mcols(x)[["regionIndex"]] <- DataFrame(lvl1=c(1,1,2,2,3,3,4,4,5,5))
	mcols(x)[["regionIndexName"]] <- DataFrame(lvl1=c(1,1,2,2,3,3,4,4,5,5))

	#prepeare expected data
	exp <- granges(a)[c(1,3,2,3,3)]
	end(exp) <-	start(exp)[c(3,3,3,3,1)]
	names(exp) <- 1:5
	mcols(exp) <- NULL
	mcols(exp)[["ASEsetIndex"]] <- IntegerList(list(c(1,2),c(3,2),c(2,2),c(3,2),c(1,3)))
	mcols(exp)[["regionIndex"]] <- DataFrame(lvl1=1:5)
	mcols(exp)[["regionIndexName"]] <- DataFrame(lvl1=as.character(c(1,2,3,4,5)))

	#run tests
	res <- .makeRegionGRangesFromASEsetWithRegionIndex(x)
	
	#test equality
    expect_that(exp, equals(res))

})
