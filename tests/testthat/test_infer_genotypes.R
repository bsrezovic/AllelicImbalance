
context("inferGenotypes ASEset")

test_that("correct inference of genotypes from ASEset return.class='matrix'", {

	#prepare testdata
	data(ASEset)	
	x <- ASEset[1:3,c(1:9)]
		
	res <- inferGenotypes(x,strand="*", return.class = "matrix",
		threshold.frequency = 0.2,
		threshold.count.sample = 0
	)

    expect_that(rownames(res), equals(rownames(x)))
    expect_that(colnames(res), equals(colnames(x)))
    expect_that(res,equals(matrix(
				c("A/C", "A/A", "G/G", "A/A", "A/A", "G/G", "A/A", "A/A",
			   NA, "A/A", "A/A", "G/G", "A/A", "A/A", "G/G", "A/C",
			   "A/A", "G/A", "A/C", "A/T", "G/G", "A/A", "A/A", "G/G",
			   NA, NA, NA),
				nrow=nrow(x), dimnames=list(rownames(x),colnames(x)),
				byrow=FALSE)))

})



