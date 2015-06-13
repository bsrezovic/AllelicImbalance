
#context("calculate reference fraction from ASEset")
#
#test_that(paste("calculate reference fraction from ASEset return.class='matrix'"), {
#
#	#prepare testdata
#	data(ASEset)	
#	x <- ASEset
#
#
#	#make this when the new reference fasta file is in place
#	res <- refFraction(x,strand="*", threshold.count.sample=1, inferGenotypes=TRUE) 
#
#	#remove one dim
#	res <- res[,,1]
#
#    expect_that(rownames(res), equals(rownames(x)))
#    expect_that(colnames(res), equals(colnames(x)))
#    expect_that(as.vector(res)[1:3], equals(c(0.8000000, 0.8076923, NaN)))
#    expect_that(round(as.vector(res)[57:60],7), equals(c(NaN, 0.7656250, 0.8163265, NaN)))
#
#})
#
#

