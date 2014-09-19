
context("inferGenotypes ASEset")

test_that("correct inference of genotypes from ASEset return.class='matrix'", {

	#prepare testdata
	data(ASEset)	
	x <- ASEset[1:3,c(1:9)]
		
	res <- inferGenotypes(x,strand="*", return.class = "matrix",
		threshold.frequency = 0.05,
		threshold.count.sample = 0
	)

    expect_that(rownames(res), equals(rownames(x)))
    expect_that(colnames(res), equals(colnames(x)))
    expect_that(res,equals(matrix(
				c("T/G", "C/G", "G/G", "T/T", "C/C", "G/G", "T/T", "C/C",
			   NA, "T/T", "C/C", "G/G", "T/T", "C/C", "G/G", "T/G",
			   "C/G", "G/A", "T/G", "C/G", "G/G", "T/T", "C/C", "G/G",
			   NA, NA, NA),
				nrow=nrow(x), dimnames=list(rownames(x),colnames(x)),
				byrow=FALSE)))

})



