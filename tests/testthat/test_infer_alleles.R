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

test_that("test inferAltAllele", {

	#prepare testdata
	ref <- c("G","G","C")
	ar <- matrix(c("G", "T", "C", "A", "G", "G"), ncol=2)
		
	res <- pickSecondMostExpressedAllele(ref, ar)
	exp <- c("A", "T", "G")

    expect_that(res, equals(exp))

})





