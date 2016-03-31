
context("internal functions for phase2genotype")

test_that(paste("checking .phaseArray2genotypeArray"), {

	#####################
	# Test 1
	#####################
	#A matrix with both alleles present in all samples (simplest case)
	#Dim: 3 SNPs 4 Samples
	ref <- c("A", "G", "C")
	alt <- c("T", "C", "G")
			
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
	dimnames(arr) <- list(NULL, NULL, NULL )
						  
	#shift 1 to 0 (so that ref equals 0 and alt equals 1)
	tf <- arr[,,1:2] == 1
	arr[,,1:2][tf] <- 0
	arr[,,1:2][!tf] <- 1

	#prepare expected data
	e1 <- c("A", "T", "A", "T",
			"C", "C", "G", "C",
			"G", "G", "G", "G")
	e2 <- c("T", "T", "A", "T",
			"G", "C", "C", "C",
			"C", "G", "C", "G")
	e3 <- c("|", "/", "|", "/",
			"/", "/", "|", "/",
			"|", "/", "/", "/")

	exp <- aperm(array(c(e1, e2, e3), c(4,3,3)), c(2,1,3))
	dimnames(exp) <- list(NULL, NULL, c("mat","pat", "phased"))
			
	#run tests
	res <- .phaseArray2genotypeArray(arr, ref, alt)
	
	#test equality
    expect_that(exp, equals(res))

	#####################
	# Test 2
	#####################
	#A matrix with both alleles present in all samples (NA case)
	#Dim: 3 SNPs 4 Samples
	ref <- c("A", "G", "C")
	alt <- c("T", "C", "G")
			
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
	dimnames(arr) <- list(NULL, NULL, NULL )

	#shift 1 to 0 (so that ref equals 0 and alt equals 1)
	tf <- arr[,,1:2] == 1
	arr[,,1:2][tf] <- 0
	arr[,,1:2][!tf] <- 1

	#prepare expected data
	e1 <- c("A", "T", "A", "T",
			NA, "C", "G", "C",
			"G", "G", "G", "G")
	e2 <- c("T", "T", "A", "T",
			"G", "C", "C", "C",
			"C", "G", "C", NA)
	e3 <- c("|", "/", "|", "/",
			"/", "/", "|", "/",
			"|", "/", "/", "/")

	exp <- aperm(array(c(e1, e2, e3), c(4,3,3)), c(2,1,3))
	dimnames(exp) <- list(NULL, NULL, c("mat","pat", "phased"))
			
	#run tests
	res <- .phaseArray2genotypeArray(arr, ref, alt)
	
	#test equality
    expect_that(exp, equals(res))
})

test_that(paste("checking .genotypeArray2genotypeMatrix"), {

	#####################
	# Test 1
	#####################
	#A matrix with both alleles present in all samples (simplest case)
	#Dim: 3 SNPs 4 Samples
			
	g1 <- c("A", "T", "A", "T",
			"C", "C", "G", "C",
			"G", "G", "G", "G")
	g2 <- c("T", "T", "A", "T",
			"G", "C", "C", "C",
			"C", "G", "C", "G")
	g3 <- c("|", "/", "|", "/",
			"/", "/", "|", "/",
			"|", "/", "/", "/")

	arr <- aperm(array(c(g1, g2, g3), c(4,3,3)), c(2,1,3))
	dimnames(arr) <- list(NULL, NULL, c("mat","pat", "phased"))

	#prepare expected data
	e1 <- c("A|T", "T/T", "A|A", "T/T",
			"C/G", "C/C", "G|C", "C/C",
			"G|C", "G/G", "G/C", "G/G")

	exp <- matrix(e1,ncol=4,  byrow=TRUE)
						  
	#run tests
	res <- .genotypeArray2genotypeMatrix(arr)
	
	#test equality
    expect_that(exp, equals(res))

	#####################
	# Test 2
	#####################
	#A matrix with both alleles present in all samples (NA case)
	#Dim: 3 SNPs 4 Samples
	g1 <- c("A", "T", "A", "T",
			NA, "C", "G", "C",
			"G", "G", "G", "G")
	g2 <- c("T", "T", "A", "T",
			"G", "C", "C", "C",
			"C", "G", "C", NA)
	g3 <- c("|", "/", "|", "/",
			"/", "/", "|", "/",
			"|", "/", NA, "/")

	arr <- aperm(array(c(g1, g2, g3), c(4,3,3)), c(2,1,3))
	dimnames(arr) <- list(NULL, NULL, c("mat","pat", "phased"))

	#prepare expected data
	e1 <- c("A|T", "T/T", "A|A", "T/T",
			NA, "C/C", "G|C", "C/C",
			"G|C", "G/G", "G/C", NA)

	exp <- matrix(e1,ncol=4,  byrow=TRUE)
						  
	#run tests
	res <- .genotypeArray2genotypeMatrix(arr)
	
	#test equality
    expect_that(exp, equals(res))

})


