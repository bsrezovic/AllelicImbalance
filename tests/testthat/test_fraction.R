
context("test fraction functionality")

test_that(paste("test fraction default"), {

	#prepare testdata
	data(ASEset)	
	x <- ASEset

	res <- t(fraction(x))
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
	res <- t(fraction(a, top.fraction.criteria="phase"))

	#check 
    expect_that(colnames(res), equals(rownames(a)))
    expect_that(rownames(res), equals(colnames(a)))
    expect_that(as.vector(res)[1:5], equals(c(0.2, 1.0, 1.0, 0.0, 1.0)))
    expect_that(as.vector(res)[57:60], equals(c(NaN,0,0,0)))

})

test_that(paste("checking .returnMaternalPhaseFrequency"), {

	#####################
	# Test 1 - 
	#####################
	
	#add one region that is missing
	ar <- aperm(array(c(T,T,F,F,F,F,F,F,F,F,F,T), dim=c(3,4,2)), c(1,3,2))
	fr <- array(c(0.3,0.2,0.3,0.2,0.8,0.2,0.3,0.7,0.4,0.2,0.1,0.3), dim=c(3,2,4))
	
	#prepare phase matrix
	ph <- array(c(1,1,0,0,1,0),dim=c(3,2))
	#prepare reference allele matrix
	rf <- array(c(0,0.2,0.8,0,0,1),dim=c(3,2))

	#prepeare expected data
	exp <- array(c(1.0,0.8,0.8,0,1,1), dim=c(3,2))

	#run tests
	res <- .returnMaternalPhaseFrequency(ph,rf)
	
	#test equality
    expect_equal(exp, res)


})


