context("positions from CIGAR")

test_that("positions including deletions are calculated correctly ", {

	#position on the ref-genome
	refPOS <- c(
		5,4,6 #deletions
	)
	#the actual position in the read
	readPOS <- c(
		4,0,0 #deletions
	)
	CIGARS <- c(
	"3M1D6M","3M1D1M1D5M","3M1D1M1D5M" #deletions
	)
	cl <- cigarToRleList(CIGARS)
	funPOS <- mapply(realCigarPosition,cl,refPOS)
	#check if equals
	
	expect_that(funPOS, equals(readPOS))
})

test_that("positions including insertions are calculated correctly ", {

	#position on the ref-genome
	refPOS <- c(
		5,4,6, #insertions
		3
	)
	#the actual position in the read
	readPOS <- c(
		8,5,8, #insertions
		6
	)
	CIGARS <- c(
		"3M3I6M","3M1I1M2I5M","3M1I1M1I5M", #insertion
		"2M3I2M2I4M"
	)
	cl <- cigarToRleList(CIGARS)
	funPOS <- mapply(realCigarPosition,cl,refPOS)

	#check if equals
	expect_that(funPOS, equals(readPOS))

})


test_that("positions including introns are calculated correctly ", {
	#position on the ref-genome
	refPOS <- c(
		20,204,98 #

	)
	#the actual position in the read
	readPOS <- c(
		-1,4,-1 #
	)
	CIGARS <- c(
	"3M200N97M","3M200N97M","97M200N3M" #introns
	)

	cl <- cigarToRleList(CIGARS)
	funPOS <- mapply(realCigarPosition,cl,refPOS)

	#check if equals
	expect_that(funPOS, equals(readPOS))


})

