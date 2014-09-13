context("inferAlleles ASEset")

test_that("correct inference of alleles from ASEset return.type='all'", {

	#prepare testdata
	data(ASEset)	
	x <- ASEset[1:3,1:3]
		
	res <- inferAlleles(x, return.type = "all",
		threshold.frequency = 0.05, inferOver="allSamples" )

    expect_that(rownames(res), equals(rownames(x)))
    expect_that(colnames(res), equals(c("uni","bi","tri","quad")))
    expect_that(matrix(res,ncol=4), 
		equals(matrix(c(rep(TRUE,5),rep(FALSE,7)),ncol=4,byrow=FALSE)))



})



