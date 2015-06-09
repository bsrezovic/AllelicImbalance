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


