
context("filter heterozygots from ASEset")

test_that(paste("filter heterozygots from ASEset"), {

	#prepare testdata
	data(ASEset)	
	x <- ASEset
	genotype(x) <- inferGenotypes(x)

	res <- hetFilt(x)

    expect_that(rownames(res), equals(rownames(x)))
    expect_that(colnames(res), equals(colnames(x)))
    expect_that(as.vector(res)[1:12], equals(c(
		TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, 
		FALSE, FALSE, NA, FALSE, FALSE, FALSE)))

})



