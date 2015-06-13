
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

