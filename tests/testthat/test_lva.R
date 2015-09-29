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

test_that(paste("checking .lvaRegressionReturnCommonParamMatrix"), {

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
	ar <- aperm(ar,c(3,1,2))
	
	grp <- t(matrix(c(1,1,3,2,2,3),nrow=2))

	#prepare expected data (values are from an execution of the linear model)
	exp <- matrix(
			c(0.23333333, 0.05, 0.31180478, 0.14433757, 0.7483315, 0.3464102, 0.5910148, 0.7877044,
			  0.03333333, 0.15, 0.06236096, 0.02886751, 0.5345225, 5.1961524, 0.6874944, 0.1210377),
			  ncol=8, nrow=2,byrow=TRUE)
	colnames(exp) <- c("est1","est2","stderr1","stderr2","tvalue1","tvalue2","pvalue1","pvalue2")

	#run tests
	res <- .lvaRegressionReturnCommonParamMatrix(ar, grp, 3)
	
	#test equality
    expect_equal(exp, res, tolerance=1e-6)

	#####################
	# Test 2 -what if there are only NAs. or only one value is non.NA
	#####################
	#
	#

	#introduce NAs in a complete row
	ar1 <- ar
	ar1[2,,] <- NA

	#introduce NAs in a complete row, but having one non-na value
	ar2 <- ar1
	ar2[2,1,] <- 1
	
	#prepare expected data (values are from an execution of the linear model)
	exp1 <- matrix(
			c(0.23333333, 0.05, 0.31180478, 0.14433757, 0.7483315, 0.3464102, 0.5910148, 0.7877044,
			  NA, NA, NA, NA, NA, NA, NA, NA),
			  ncol=8, nrow=2,byrow=TRUE)
	exp2 <- matrix(
			c(0.23333333, 0.05, 0.31180478, 0.14433757, 0.7483315, 0.3464102, 0.5910148, 0.7877044,
			  1.00000000, NaN, NaN, NaN, NA, NA, NA, NA),
			  ncol=8, nrow=2,byrow=TRUE)
	colnames(exp1) <- c("est1","est2","stderr1","stderr2","tvalue1","tvalue2","pvalue1","pvalue2")
	colnames(exp2) <- c("est1","est2","stderr1","stderr2","tvalue1","tvalue2","pvalue1","pvalue2")

	#run tests
	res1 <- .lvaRegressionReturnCommonParamMatrix(ar1, grp, 3)
	res2 <- .lvaRegressionReturnCommonParamMatrix(ar2, grp, 3)
	
	#test equality
    expect_equal(exp1, res1, tolerance=1e-6)
    expect_equal(exp2, res2, tolerance=1e-6)
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
    expect_that(exp, equals(res))

	#####################
	# Test 2 - test that same rows give same values in the result
	#####################
	# rows: SNPs
	# cols: samples
			
	mat <- c(0, 1, 1, 1,
			 1, 0, 1, 1,
			 0, 1, 1, 1,
			 0, 1, 1, 1,
			 0, 1, 1, 1)

	pat <- c(0, 1, 0, 1,
			 0, 1, 1, 0,
			 1, 1, 0, 0,
			 1, 1, 0, 0,
			 1, 1, 0, 0)

	ar <- aperm(array(c(mat, pat), c(4, 5, 2)),c(2,1,3))

	#prepare expected data (homozygotes should have 2)
	e1  <- c(2, 2, 1, 2,
			 1, 3, 2, 1,
			 3, 2, 1, 1,
			 3, 2, 1, 1,
			 3, 2, 1, 1)
	exp <- matrix(e1,ncol=4, byrow=TRUE)

	#run tests
	res <- .groupBasedOnPhaseAndAlleleCombination(ar)
	
	#test equality
    expect_that(exp, equals(res))

})

test_that(paste("checking .groupBasedOnPhaseAndAlleleCombination"), {


	#####################
	# Test 2 - test that same rows give same values in the result
	#####################
	# rows: SNPs
	# cols: samples
			
	mat <- c(0, 1, 1, 1,
			 1, 0, 1, 1,
			 0, 1, 1, 1,
			 0, 1, 1, 1,
			 0, 1, 1, 1)

	pat <- c(0, 1, 0, 1,
			 0, 1, 1, 0,
			 1, 1, 0, 0,
			 1, 1, 0, 0,
			 1, 1, 0, 0)

	ar <- aperm(array(c(mat, pat), c(4, 5, 2)),c(2,1,3))

	#prepare expected data (homozygotes should have 2)
	e1  <- c(2, 2, 1, 2,
			 1, 3, 2, 1,
			 3, 2, 1, 1,
			 3, 2, 1, 1,
			 3, 2, 1, 1)
	exp <- matrix(e1,ncol=4, byrow=TRUE)

	#run tests
	res <- .groupBasedOnPhaseAndAlleleCombination(ar)
	
	#test equality
    expect_that(exp, equals(res))

})

