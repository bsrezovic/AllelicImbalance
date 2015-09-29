context("test utils functions")

test_that(paste("test trailing slash remover"), {


	dir1 <- system.file('extdata/ERP000101_subset', package='AllelicImbalance')
	dir2 <- paste(dir1,"/",sep="")	

	files1 <- c("ERR009113.bam","ERR009147.bam", "ERR009159.bam" )
	files2 <- paste(files1,"/",sep="")	
	files3 <- paste("/", files1,sep="")	


	##prepare result string (at the moment only system specific)
	#res <- c(
	#		 paste("/mnt/kelewan/pappewaio/Documents/PHD/repos/AllelicImbalance",
	#			 "/bioCgit/AllelicImbalance/inst/extdata/ERP000101_subset/",
	#			 "ERR009113.bam",
	#			 sep=""),
	#		 paste("/mnt/kelewan/pappewaio/Documents/PHD/repos/AllelicImbalance",
	#			 "/bioCgit/AllelicImbalance/inst/extdata/ERP000101_subset/",
	#			 "ERR009147.bam",
	#			 sep=""),
	#		 paste("/mnt/kelewan/pappewaio/Documents/PHD/repos/AllelicImbalance",
	#			 "/bioCgit/AllelicImbalance/inst/extdata/ERP000101_subset/",
	#			 "ERR009159.bam",
	#			 sep="")
	#	)

	##check single file
    #expect_that(res[1], equals(.mergeDirAndFilename(dir1, files1[1])))
    #expect_that(res[1], equals(.mergeDirAndFilename(dir1, files2[1])))
    #expect_that(res[1], equals(.mergeDirAndFilename(dir1, files3[1])))

    #expect_that(res[1], equals(.mergeDirAndFilename(dir2, files1[1])))
    #expect_that(res[1], equals(.mergeDirAndFilename(dir2, files2[1])))
    #expect_that(res[1], equals(.mergeDirAndFilename(dir2, files3[1])))

	##check multiple files
    #expect_that(res, equals(.mergeDirAndFilename(dir1, files1)))
    #expect_that(res, equals(.mergeDirAndFilename(dir1, files2)))
    #expect_that(res, equals(.mergeDirAndFilename(dir1, files3)))

    #expect_that(res, equals(.mergeDirAndFilename(dir2, files1)))
    #expect_that(res, equals(.mergeDirAndFilename(dir2, files2)))
    #expect_that(res, equals(.mergeDirAndFilename(dir2, files3)))

})

test_that(paste("test lm simplifier"), {

	#######################
	# Test 1
	#######################
	y <- c(1,2,4,3,2,4,3,2)
	x <- c(1,1,1,1,2,2,2,2)
	en	<- lm(x~y)
	lst <- list(en,en,en,en)

	exp <- matrix(c(1.333333, 0.06349206, 0.5727498, 0.2040984,
				  2.327951, 0.3110855, 0.05880542, 0.7662601),
				  ncol=8, nrow=4, byrow=TRUE)
	colnames(exp) <- c("est1","est2","stderr1","stderr2","tvalue1","tvalue2","pvalue1","pvalue2")
	
	res <- .matrixFromLmListCommonParam(lst)

    expect_equal(res, exp, tolerance=1e-5)

})

test_that(paste("checking .IRangesFromIntegerList"), {

	#####################
	# Test 1 
	#####################
	
	#add one region that is missing
	idx <- IntegerList(list(c(1,2),c(3,3,4,1),c(4,4,3)))

	#prepeare expected data
	exp <- IRanges(c(1,3,7), c(2,6,9))

	#run tests
	res <- .IRangesFromIntegerList(idx)
	
	#test equality
    expect_that(exp, equals(res))

})
