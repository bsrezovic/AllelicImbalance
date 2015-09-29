context("test units for method lva and lva.internal")

test_that(paste("checking .lvaRegressionPvalue"), {

	#####################
	# Test 1
	#####################
	#
	#
	hets <- c(3, 4, 2)
	homs <- c(1, 2, 2)
	mean.fr <- c(0.2, 0.3, 0.5)
	sd.fr <- c(NA, NA, 0.3)
	mean.delta <- c(NA, NA, 0.2)
	sd.delta <- c(NA, NA, 0.2)
	ai.up <- c(1, 2, 2)
	ai.down <- c(0, 1, 2)
						  
	ar <- array(c(hets,homs,mean.fr, sd.fr, mean.delta, sd.delta, ai.up, ai.down), dim=c(3,8,2))
	grp <- matrix(c(1,1,3,2,2,3),nrow=2)

	#prepare expected data (values are from an execution of the linear model)
	exp <- list(0.7877044, 0.1210377)

	#run tests
	res <- .lvaRegressionPvalue(ar, grp, 3)
	
	#test equality
    expect_equal(exp, res, tolerance=1e-6)

})

test_that(paste("checking .groupBasedOnPhaseAndAlleleCombination"), {

	#####################
	# Test 1
	#####################
	# rows: SNPs
	# cols: samples
			
	mat <- c(0, 1, 1, 1,
			 1, 0, 1, 1,
			 0, 1, 1, 1)

	pat <- c(0, 1, 0, 1,
			 0, 1, 1, 0,
			 1, 1, 0, 0)

	ar <- aperm(array(c(mat, pat), c(4, 3, 2)),c(2,1,3))

	#prepare expected data (homozygotes should have 2)
	e1  <- c(2, 2, 1, 2,
			 1, 3, 2, 1,
			 3, 2, 1, 1)
	exp <- matrix(e1,ncol=4, byrow=TRUE)

	#run tests
	res <- .groupBasedOnPhaseAndAlleleCombination(ar)
	
	#test equality
    expect_equal(as.numeric(exp), as.numeric(res))

})


<<<<<<< HEAD

=======
>>>>>>> chopped out linkage and summary methods R files and chopped out a first batch of units and tested them
