
context("internal functions for genotype2phase")

test_that(paste("checking .splitGenotypeMatrix"), {

	#####################
	# Test 1
	#####################
	#A matrix with both alleles present in all samples (simplest case)
	#Dim: 3 SNPs 4 Samples
	a1 <- c("A", "A", "T", "T",
			"G", "G", "C", "C",
			"C", "G", "C", "G")
	a2 <- c("T", "A", "A", "T",
			"C", "C", "C", "G",
			"G", "G", "C", "C")
	mat <- matrix(paste(a1,"/",a2,sep=""), nrow=3, ncol=4, byrow=TRUE)
	exp <- aperm(array(c(a1, a2), c(4,3,2)),c(2,1,3))
			
	#run tests
	res <- .splitGenotypeMatrix(mat)
	
	#test equality
    expect_that(exp, equals(res))

	#####################
	# Test 2
	#####################
	#A matrix with two genotypes missing
	#Dim: 3 SNPs 4 Samples
	a1 <- c("A", "A", "T", "T",
			"G", "G", "C", "C",
			"C", "G", "C", "G")
	a2 <- c("T", "A", "A", "T",
			"C", "C", "C", "G",
			"G", "G", "C", "C")
	mat <- matrix(paste(a1,"/",a2,sep=""), nrow=3, ncol=4, byrow=TRUE)
	mat[c(1,6)] <- NA
	exp <- aperm(array(c(a1, a2), c(4,3,2)),c(2,1,3))
	exp[c(1,6,1+12,6+12)] <- NA
			
	#run tests
	res <- .splitGenotypeMatrix(mat)
	
	#test equality
    expect_that(exp, equals(res))

})

test_that(paste("checking .splitGenotypeCount"), {

	#####################
	# Test 1
	#####################
	#A matrix with both alleles present in all samples (simplest case)
	#Dim: 3 SNPs 4 Samples
	a1 <- c("A", "A", "T", "T",
			"G", "G", "C", "C",
			"C", "G", "C", "G")
	a2 <- c("T", "A", "A", "T",
			"C", "C", "C", "G",
			"G", "G", "C", "C")
	ar <- aperm(array(c(a1, a2), c(4,3,2)),c(2,1,3))
	exp <- matrix(c(4, 0, 0, 4, 
					0, 5, 3, 0,
					0, 4, 4, 0), nrow=3, byrow=TRUE)
	colnames(exp) <- c("A","C","G","T")
			
	#run tests
	res <- .splitGenotypeCount(ar)
	
	#test equality
    expect_that(exp, equals(res))

})

test_that(paste("checking .splitGenotypeRank"), {

	#####################
	# Test1
	#####################
	#A matrix with both alleles present in all samples (simplest case)
	#Dim: 3 SNPs 4 Samples
	mat <- matrix(c(4, 0, 0, 4, 
					0, 5, 3, 0,
					0, 4, 4, 0), nrow=3, byrow=TRUE)
	colnames(mat) <- c("A","C","G","T")
	exp <- matrix(c("A", "T", "C", "G",
					"C", "G", "A", "T",
					"C", "G", "A", "T"), nrow=3, byrow=TRUE)
	colnames(exp) <- c("r1","r2","r3","r4")
			
	#run tests
	res <- .splitGenotypeRank(mat)
	
	#test equality
    expect_that(exp, equals(res))

})



