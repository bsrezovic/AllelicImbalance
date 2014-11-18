context("count alleles from Bam File")

#test_that(paste("check that different import methods give same result"), {
#
#	#prepare testdata
#	data(GRvariants)
#	gr <- GRvariants
#	pathToDir <- system.file('inst/extdata/ERP000101_subset', package='AllelicImbalance')
#	#pathToDir <- system.file('inst/extdata/ERP000101_subset/tmp', package='AllelicImbalance')
#
#	#import variant 1
#	bam <- impBamGAL(pathToDir, gr, reduce=TRUE, readLength=36,
#					 scanBamFlag=scanBamFlag(), verbose=FALSE)
#	bam2 <- GAlignmentsList(lapply(bam,function(x){x[mcols(x)[["flag"]]%in%c(99,147,163,83)]}))
#	seqlevels(bam2) <- seqlevels(gr)
#
#	#change strand to match secondstrand
#	bam2 <- GAlignmentsList(lapply(bam2,function(x){
#					cat(names(x))
#					strand(x)[mcols(x)[["flag"]]%in%c(83,163)] <- "-"
#					strand(x)[mcols(x)[["flag"]]%in%c(99,147)] <- "+"
#					x
#				  }))
#
#	plus <- getAlleleCounts(bam2, gr, strand="+",return.class="array",verbose=F )
#	minus <- getAlleleCounts(bam2, gr, strand="-",return.class="array",verbose=F)
#
#	#import variant 2
#
#	#minus strand
#	arm1 <- countAllelesFromBam(gr, pathToDir, flag=83, verbose=FALSE)
#	arm2 <- countAllelesFromBam(gr, pathToDir, flag=163, verbose=FALSE)
#	arm <- arm1 + arm2
#
#	#plus strand
#	arp1 <- countAllelesFromBam(gr, pathToDir, flag=99, verbose=FALSE)
#	arp2 <- countAllelesFromBam(gr, pathToDir, flag=147, verbose=FALSE)
#	arp <- arp1 + arp2
#
#	#test equality
#    expect_that(rownames(plus), equals(rownames(arp)))
#    expect_that(colnames(plus), equals(colnames(arp)))
#    expect_that(rownames(minus), equals(rownames(arm)))
#    expect_that(colnames(minus), equals(colnames(arm)))
#	expect_identical(plus, arp)
#	expect_identical(minus, arm)
#})



