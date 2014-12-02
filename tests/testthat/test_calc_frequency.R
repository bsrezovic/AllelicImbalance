
context("frequency ASEset")

test_that("frequency list from ASEset is calculated correctly", {
	
	#prepare testdata
	data(ASEset)	
	x <- ASEset[1:2,1:3]
	
	#test function
	fr <- frequency(x, return.class="list", threshold.count.sample = 1)

    # check 
    expect_that(names(fr), equals(rownames(x)))
    expect_that(fr[[1]], 
		equals(matrix(c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.0, 0.0, 0.8, 1.0, 1.0),
		ncol=4, dimnames=list(colnames(x), x@variants))))
    expect_that(as.numeric(fr[[2]][,2]), equals(c(0.8076923, 1.0000000, 1.0000000 )))
    expect_that(round(as.numeric(fr[[2]][,3]),7), equals(c(0.1923077, 0.0000000, 0.0000000 )))
    expect_that(as.numeric(apply(fr[[1]],1,sum)), equals(c(1,1,1)))
    expect_that(as.numeric(apply(fr[[2]],1,sum)), equals(c(1,1,1)))

})

test_that("frequency array from ASEset is calculated correctly", {
	
	#prepare testdata
	data(ASEset)	
	x <- ASEset[1:3,1:3]
	
	#test function
	fr <- frequency(x, return.class="array", threshold.count.sample = 1)

    # check 
    expect_that(dimnames(fr)[[1]], equals(rownames(x)))
    expect_that(dimnames(fr)[[2]], equals(colnames(x)))
    expect_that(dimnames(fr)[[3]], equals(x@variants))
    expect_that(fr[3,3,3], equals(NaN))
    expect_that(as.numeric(apply(fr[,1,],1,sum)), equals(c(1,1,1)))
    expect_that(as.numeric(apply(fr[,2,],1,sum)), equals(c(1,1,1)))
    expect_that(as.numeric(apply(fr[,3,],1,sum)), equals(c(1,1,NaN)))

})


