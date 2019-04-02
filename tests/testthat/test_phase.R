
context("test phasing functionality")

test_that(paste("test store and access of phase"), {

	#prepare testdata
	data(ASEset)	
	x <- ASEset

	set.seed(1)
	p1 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(x)*ncol(x)),nrow=nrow(x), ncol(x))
	p2 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(x)*ncol(x)),nrow=nrow(x), ncol(x))
	p <- matrix(paste(p1,sample(c("|","|","/"), size=nrow(x)*ncol(x), replace=TRUE), p2, sep=""),
		nrow=nrow(x), ncol(x))

	#store matrix
	phase(x) <- p

	#store array
	phase(x) <- phaseMatrix2Array(p, dimnames=dimnames(x))

	#access matrix
	res <- phase(x)

	#access array
	res2 <- phase(x, return.class="array")

	#check 
    expect_that(rownames(res), equals(rownames(x)))
    expect_that(colnames(res), equals(colnames(x)))
    expect_that(as.vector(res)[1:4], equals(c("1|1","0/0","1/0","1/0")))
    expect_that(as.vector(res)[57:60], equals(c("0/0","0|1","0|0","1/0")))

    expect_that(rownames(res2), equals(rownames(x)))
    expect_that(colnames(res2), equals(colnames(x)))
    expect_that(as.vector(res2)[1:4], equals(c(1,0,1,1)))
    expect_that(as.vector(res2)[177:180], equals(c(0,1,1,0)))

})

