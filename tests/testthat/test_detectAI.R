context("detect AI from ASEset")

test_that(paste("detectAI from ASEset results and accessors, return.class='DetectedAI',", 
				"multiple samples, multipliple SNPs"), {

	#prepare testdata
	data(ASEset)	
	x <- ASEset

	#make this when the new reference fasta file is in place
	res <- detectAI(x,strand="*", threshold.count.sample=1) 

	#check row and colnames
    expect_that(rownames(res), equals(rownames(x)))
    expect_that(colnames(res), equals(colnames(x)))

	#check strand
    expect_that(rownames(res), equals(rownames(x)))

	#check reference.frequency
    expect_that(as.vector(referenceFrequency(res))[1:3], equals(c(0.8000000, 0.8076923, NaN)))
    expect_that(round(referenceFrequency(res)[57:60],7), equals(c(NaN, 0.7656250, 0.8163265, NaN)))

	#check threshold.frequency
    expect_that(thresholdFrequency(res)[1:3], equals(c(TRUE,TRUE,NA)))
    expect_that(thresholdFrequency(res)[57:60], equals(c(NA, TRUE, TRUE, NA)))
	
	#check threshold.count.sample
    expect_that(thresholdCountSample(res)[1:3], equals(c(TRUE,TRUE,FALSE)))
    expect_that(thresholdCountSample(res)[57:60], equals(c(FALSE, TRUE, TRUE, FALSE)))

	#check threshold.delta.frequency
    expect_that(thresholdDeltaFrequency(res)[1:3], equals(c(TRUE,TRUE,NA)))
    expect_that(thresholdDeltaFrequency(res)[57:60], equals(c(NA, TRUE, TRUE, NA)))

	#check threshold.pvalue
    expect_that(thresholdPvalue(res)[1:3], equals(c(TRUE,TRUE,NA)))
    expect_that(thresholdPvalue(res)[57:60], equals(c(NA, TRUE, TRUE, NA)))

})



