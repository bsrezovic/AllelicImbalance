
context("test fraction functionality")

test_that(paste("test fraction default"), {

	#prepare testdata
	data(ASEset)	
	x <- ASEset

	res <- fraction(x)

	#check 
    expect_that(colnames(res), equals(rownames(x)))
    expect_that(rownames(res), equals(colnames(x)))
    expect_that(as.vector(res)[1:3], equals(c(0.8, 1.0, 1.0)))
    expect_that(as.vector(res)[57:60], equals(c(NaN,1,1,1)))

})


test_that(paste("test fraction usePhase=TRUE"), {

	#prepare testdata
	#load data
	data(ASEset) 
    a <- ASEset
  	                                                                                              
	set.seed(1)
    #example phase matrix 
    p1 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
    p2 <- matrix(sample(c(1,0),replace=TRUE, size=nrow(a)*ncol(a)),nrow=nrow(a), ncol(a))
    pha <- matrix(paste(p1,sample(c("|","|","/"), size=nrow(a)*ncol(a), replace=TRUE), p2, sep=""),
    	nrow=nrow(a), ncol(a))
  	                                                                                              
	#infer alt allele
	mcols(a)[["alt"]] <- inferAltAllele(a)

	#store
	phase(a) <- pha
	
	#check 
	res <- fraction(a, top.fraction.criteria="phase")

	#check 
    expect_that(colnames(res), equals(rownames(a)))
    expect_that(rownames(res), equals(colnames(a)))
    expect_that(as.vector(res)[1:5], equals(c(0.2, 1.0, 1.0, 0.0, 1.0)))
    expect_that(as.vector(res)[57:60], equals(c(NaN,0,0,0)))

})

