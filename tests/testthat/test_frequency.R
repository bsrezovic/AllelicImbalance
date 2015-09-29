context("test units for method frequency")

test_that(paste("checking .calcFrequencyFromAlleleCounts"), {

	#####################
	# Test 1
	#####################
	#Dim: 3 SNPs 4 Samples
			
	a1 <- c(1, 2, 1, 2,
			2, 3, 2, 1,
			3, 4, 1, 2)
	a2 <- c(0, 2, 1, 3,
			2, 3, 2, 4,
			1, 0, 1, 6)
	a3 <- c(2, 2, 1, 0,
			4, 2, 2, 0,
			0, 1, 0, 0)
	a4 <- c(1, 2, 1, 0,
			2, 0, 4, 0,
			1, 0, 0, 0)

	arr <- aperm(array(c(a1, a2, a3, a4), dim=c(4,3,4)), c(2,1,3))
	dimnames(arr) <- list(paste("SNP",1:3,sep=""), paste("sample",1:4,sep=""), paste("allele",1:4,sep=""))
						  
	#prepare expected data
	e1 <- c(0.25, 0.250, 0.25, 0.40,
			0.20, 0.375, 0.20, 0.20,
			0.60, 0.800, 0.50, 0.25)
	e2 <- c(0.0, 0.250, 0.25, 0.60,
			0.2, 0.375, 0.20, 0.80,
			0.2, 0.000, 0.50, 0.75)
	e3 <- c(0.5, 0.25, 0.25, 0,
			0.4, 0.25, 0.20, 0,
			0.0, 0.20, 0.00, 0)
	e4 <- c(0.25, 0.25, 0.25, 0,
			0.20, 0.00, 0.40, 0,
			0.20, 0.00, 0.00, 0)

	exp <- aperm(array(c(e1, e2, e3, e4), dim=c(4,3,4)), c(2,1,3))
	dimnames(exp) <- list(paste("SNP",1:3,sep=""), paste("sample",1:4,sep=""), paste("allele",1:4,sep=""))
						  

	#run tests
	res <- .calcFrequencyFromAlleleCounts(arr, 1)
	
	#test equality
    expect_that(exp, equals(res))

})



