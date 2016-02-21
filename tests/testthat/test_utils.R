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

test_that(paste("checking .arrayFromAlleleVector"), {

	#####################
	# Test 1 
	#####################
	
	#add one region that is missing
	var <- c("A", "T", "G", "C")
	sel <- c("A", "A","C")
	nc <- 2

	#prepeare expected data
	exp <- array(c(c(T,T,F,T,T,F,F,F,F,F,F,F),
				   c(F,F,F,F,F,F,F,F,T,F,F,T))
				   , dim=c(3,2,4))

	#run tests
	res <- .arrayFromAlleleVector(var, sel, nc)
	
	#test equality
    expect_equal(exp, res)


	#####################
	# Test 2 - something the unit test 1, did not catch
	#####################
	
	#add one region that is missing
	var <- c("A", "C", "G", "T")
	sel <- c("T", "C","G")
	nc <- 20

	#prepeare expected data
	ma <- array(c(F,F,F,T,F,T,F,F,F,F,T,F), dim=c(4,3))
	exp <- aperm(array(ma,dim=c(nrow(ma),ncol(ma),nc)),c(2,3,1) )
 
	#run tests
	res <- .arrayFromAlleleVector(var, sel, nc)
	
	#test equality
    expect_equal(exp, res)


})

test_that(paste("checking .subsetArrayToMatrix"), {

	#####################
	# Test 1 - when array to be subsetted is the frequency array
	#####################
	
	#add one region that is missing
	ar2 <- aperm(array(c(T,T,F,F,F,F,F,F,F,F,F,T), dim=c(3,4,2)), c(1,3,2))
	ar <- array(c(0.3,0.2,0.3,0.2,0.8,0.2,0.3,0.7,0.4,0.2,0.1,0.3), dim=c(3,2,4))

	#prepeare expected data
	exp <- array(c(0.3,0.2,0.4,0.2,0.8,0.3), dim=c(3,2))

	#run tests
	res <- .subsetArrayToMatrix(ar, ar2)
	
	#test equality
    expect_equal(exp, res)

	#####################
	# Test 2 - when array to be subsetted is the count array
	#####################
	
	#add one region that is missing
	ar2 <- aperm(array(c(T,T,F,F,F,F,F,F,F,F,F,T), dim=c(3,4,2)), c(1,3,2))
	ar <- array(c(10,20,30,40,50,60,11,22,33,44,55,66), dim=c(3,2,4))

	#prepeare expected data
	exp <- array(c(10,20,33,40,50,66), dim=c(3,2))

	#run tests
	res <- .subsetArrayToMatrix(ar, ar2)
	
	#test equality
    expect_equal(exp, res)

})

test_that(paste("checking .expandMatrixToArray"), {

	#####################
	# Test 1 - 
	#####################
	
	#prepeare expected data
	mat <- array(c(0.3,0.2,0.4,0.2,0.8,0.3), dim=c(3,2))
	len <- 4

	#add one region that is missing
	exp <- array(mat, dim=c(3,2,4))

	#run tests
	res <- .expandMatrixToArray(mat, len)
	
	#test equality
    expect_equal(exp, res)

})

test_that(paste("checking .verboseCoerceToCharacter"), {

	#####################
	# Test 1 - 
	#####################
	
	#prepeare expected data
	vec <- c("A","T","G")

	#add one region that is missing
	exp <- c("A","T","G")

	#run tests
	res <- .verboseCoerceToCharacter(vec)
	
	#test equality
    expect_equal(exp, res)

	#####################
	# Test 1 - test a class that is not character
	#####################

	vec <- c(1,2,3)
	#add one region that is missing
	exp <- c("1","2","3")

	#test 
    expect_warning(res <- .verboseCoerceToCharacter(vec))
    expect_equal(exp, res)

})

test_that(paste("checking .Na2False"), {

	#####################
	# Test 1 - vector
	#####################
	
	#prepeare expected data
	vec <- c("A",NA,"G")
	#add one region that is missing
	exp <- c("A",FALSE,"G")
	#run tests
	res <- .Na2False(vec)
	#test equality
    expect_equal(exp, res)

	#####################
	# Test 2 - test a matrix with numbers
	# numbers will not give the expected output
	#####################

	vec <- matrix(c(1,2,3,NA,5,6), ncol=2)
    expect_error(.Na2False(vec))

	#####################
	# Test 2 - test a matrix with numbers
	# numbers will not give the expected output
	#####################
	vec <- matrix(c(T,F,T,NA,F,T), ncol=2)
	exp <- matrix(c(T,F,T,F,F,T), ncol=2)
	res <- .Na2False(vec)
	expect_equal(exp, res)

})

