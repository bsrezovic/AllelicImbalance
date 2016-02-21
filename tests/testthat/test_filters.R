
context("internal functions for hetFilt")

test_that(paste("checking .heterozygozityFromPhaseArray"), {

	#####################
	# Test 1
	#####################
	#A matrix with both alleles present in all samples (simplest case)
	#Dim: 3 SNPs 4 Samples
	p1 <- c(1, 0, 1, 0,
			0, 0, 1, 0,
			0, 0, 0, 0)
	p2 <- c(0, 0, 1, 0,
			1, 0, 0, 0,
			1, 0, 1, 0)
	p3 <- c(1, 0, 1, 0,
			0, 0, 1, 0,
			1, 0, 0, 0)

	arr <- aperm(array(c(p1, p2, p3), dim=c(4,3,3)), c(2,1,3))
						  
	#prepare expected data
	e1 <- c(TRUE, FALSE, FALSE, FALSE,
			TRUE, FALSE, TRUE, FALSE,
			TRUE, FALSE, TRUE, FALSE)

	exp <- matrix(e1,ncol=4,  byrow=TRUE)
			
	#run tests
	res <- .heterozygozityFromPhaseArray(arr)
    expect_that(exp, equals(res))

	#####################
	# Test 2
	#####################
	#A matrix with both alleles present in all samples (NA case)
	#Dim: 3 SNPs 4 Samples
			
	p1 <- c(1, 0, 1, 0,
			NA, 0, 1, 0,
			0, 0, 0, 0)
	p2 <- c(0, 0, 1, 0,
			1, 0, 0, 0,
			1, 0, 1, NA)
	p3 <- c(1, 0, 1, 0,
			0, 0, 1, 0,
			1, NA, 0, 0)

	arr <- aperm(array(c(p1, p2, p3), dim=c(4,3,3)), c(2,1,3))
						  
	#prepare expected data
	e1 <- c(TRUE, FALSE, FALSE, FALSE,
			NA, FALSE, TRUE, FALSE,
			TRUE, FALSE, TRUE, NA)

	exp <- matrix(e1,ncol=4,  byrow=TRUE)
			
	#run tests
	res <- .heterozygozityFromPhaseArray(arr)
    expect_that(exp, equals(res))

})

context("internal functions for minCountFilt")

test_that(paste("checking .toKeepMatrixMinCountFilterEach"), {

	#####################
	# Test 1
	#####################
	

	ref <-	c("T","C","G")
	alt <- c("G","G","A")
	nc <- 7
	thr <- 10
	var <- c("A","C","G","T")

	p1 <- c(
	  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 10,  0,  0,
	  0,  0, 21,  0,  0, 9,  0,  0,  3,  0,  0,  1,  0,  0, 16,  0,  0, 56,  0,  0,
	  36,  0,  5,  5, 13,  0,  0,  3,  0,  0, 0,  0,  0,  2,  0,  0,  4, 19,  9, 13,
	  17, 20, 45, 20,  0,  0, 10,  0,  0,  3,  0,  0,  1,  0,  0, 15,  0,  0, 59,  0,
	  0, 38,  0,  0
	) 
	ac <- array(p1,dim=c(3,nc,4))
						  
	#prepare expected data
	e1 <- c(
		FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
		FALSE, FALSE, FALSE,  TRUE, FALSE, FALSE,  TRUE,  TRUE, FALSE

		)
	exp <- array(e1, dim=c(3,7))

			
	#run tests
	res <- .toKeepMatrixMinCountFilterEach(ac,ref,alt,nc,thr,var)
    expect_that(exp, equals(res))


	#####################
	# Test 2 - make sure the test throws error if input has dimension missmatch
	#####################
		
	expect_error(.toKeepMatrixMinCountFilterEach(ac[1:2,,], ref, alt, nc, thr, var))
	expect_error(.toKeepMatrixMinCountFilterEach(ac[1:2,,1,drop=FALSE], ref, alt, nc, thr, var))
	expect_error(.toKeepMatrixMinCountFilterEach(ac[,,1:3], ref, alt, nc, thr, var))
	expect_error(.toKeepMatrixMinCountFilterEach(ac[,1:4,], ref, alt, nc, thr, var))
	expect_error(.toKeepMatrixMinCountFilterEach(ac, ref[1:2], alt, nc, thr, var))
	expect_error(.toKeepMatrixMinCountFilterEach(ac, ref , alt[1:2], nc, thr, var))
	expect_error(.toKeepMatrixMinCountFilterEach(ac, ref , alt, nc[c(1,1)], thr, var))
	expect_error(.toKeepMatrixMinCountFilterEach(ac, ref , alt, nc, thr[c(1,1)], var))

})
test_that(paste("checking .toKeepMatrixMinCountFilterAll"), {

	#####################
	# Test 1
	#####################
	

	nc <- 7
	thr <- 10

	p1 <- c(
	  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 10,  0,  0,
	  0,  0, 21,  0,  0, 9,  0,  0,  3,  0,  0,  1,  0,  0, 16,  0,  0, 56,  0,  0,
	  36,  0,  5,  5, 13,  0,  0,  3,  0,  0, 0,  0,  0,  2,  0,  0,  4, 19,  9, 13,
	  17, 20, 45, 20,  0,  0, 10,  0,  0,  3,  0,  0,  1,  0,  0, 15,  0,  0, 59,  0,
	  0, 38,  0,  0
	) 
	ac <- array(p1,dim=c(3,nc,4))
						  
	#prepare expected data
	e1 <- c(
			  TRUE,  TRUE,  TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
			  TRUE,  TRUE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE
		)
	exp <- array(e1, dim=c(3,7))

	#run tests
	res <- .toKeepMatrixMinCountFilterAll(ac,thr)
    expect_that(exp, equals(res))


	#####################
	# Test 2 - make sure the test throws error if input has dimension missmatch
	#####################
		
	#expect_error(.toKeepMatrixMinCountFilterEach(ac[1:2,,], ref, alt, nc, thr, var))

})
