
context("internal functions for phaseArray2phaseMatrix")

test_that(paste("checking .mergePhaseArray2phaseMatrix"), {

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
	e1 <- c("1|0", "0/0", "1|1", "0/0",
			"0/1", "0/0", "1|0", "0/0",
			"0|1", "0/0", "0/1", "0/0")

	exp <- matrix(e1,ncol=4,  byrow=TRUE)
			
	#run tests
	res <- .mergePhaseArray2phaseMatrix(arr)
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
	e1 <- c("1|0", "0/0", "1|1", "0/0",
			NA, "0/0", "1|0", "0/0",
			"0|1", "0/0", "0/1", NA)

	exp <- matrix(e1,ncol=4,  byrow=TRUE)
			
	#run tests
	res <- .mergePhaseArray2phaseMatrix(arr)
    expect_that(exp, equals(res))
})

