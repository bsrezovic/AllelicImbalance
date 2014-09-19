#' detectAI
#' 
#' detection of AllelicImbalance 
#' 
#' threshold.frequency is the least fraction needed to classify as bi tri or
#' quad allelic SNPs. If 'all' then all of bi tri and quad allelic SNPs will use the same
#' threshold. Everything under the treshold will be regarded as noise. 'all' will return 
#' a matrix with snps as rows and uni bi tri and quad will be columns. For this function
#' Anything that will return TRUE for tri-allelicwill also return TRUE for uni and bi-allelic
#' for the same SNP an Sample. 
#'
#' return.type 'ref' return only AI when reference allele is more expressed. 'alt' return only
#' AI when alternative allele is more expressed or 'all' for both 'ref' and 'alt' alleles.
#' Reference allele is the one present in the reference genome on the forward strand.
#'
#' min.delta.frequency and function.test will use the value in mapBias(x) as expected value. 
#'
#' function.test will use the two most expressed alleles for testing. Make therefore sure there
#' are no tri-allelic SNPs or somatic mutations among the SNPs in the ASEset. 
#'
#' @name detectAI
#' @rdname detectAI
#' @aliases detectAI
#' detectAI,ASEset-method 
#' @docType methods
#' @param x ASEset
#' @param strand strand to infer from
#' @param threshold.count.sample least amount of counts to try to infer allele
#' @param threshold.frequency least fraction to classify (see details)
#' @param return.type 'ref' ,'alt' or default: 'all'
#' @param return.class class to return (atm only class 'logical')
#' @param min.delta.frequency minimum of frequency difference from 0.5 (or mapbias adjusted value)
#' @param max.pvalue pvalue over this number will be filtered out
#' @param function.test At the moment the only available option is 'binomial.test'
#' @param ... internal arguments 
#' @author Jesper R. Gadin
#' @keywords infer
#' @examples
#' 
#' #load example data
#' data(ASEset)
#' a <- ASEset
#'
#' dai <- detectAI(a)
#' 
#'
#' @rdname detectAI
#' @export
setGeneric("detectAI", function(x, ...){
    standardGeneric("detectAI")
})

#' @rdname detectAI
#' @export
setMethod("detectAI", signature(x = "ASEset"), function(x, 
	return.class = "logical", return.type="all", strand = "*",
	threshold.frequency=0.05, threshold.count.sample=5,
	min.delta.frequency=0.05, max.pvalue=0.05,
	function.test="binom.test") {

	fr <- refFraction(x, strand=strand,
			  threshold.count.sample=threshold.count.sample)
	fr[fr < threshold.frequency] <- NaN
	fr[fr > (1-threshold.frequency)] <- NaN

	#survive min delta freq
	biasmatRef <- matrix(aperm(mapBias(x, return.class="array"),c(3,2,1))[aperm(array(matrix(
		x@variants, ncol=length(x@variants),
		 nrow=nrow(x), byrow=TRUE) == mcols(x)[,"ref"]
	 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
	],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))

	fr2 <- fr
	fr2[is.na(fr2)] <- biasmatRef[is.na(fr2)] 

	if(any(fr2<biasmatRef)){
		fr3 <- fr2
		fr3[fr2<biasmatRef] <- biasmatRef[fr2<biasmatRef] - fr[fr2<biasmatRef]
		fr[fr3 < min.delta.frequency] <- NaN
	}

	if(any(fr2>biasmatRef)){
		fr3 <- fr2
		fr3[fr2>biasmatRef] <- fr[fr2>biasmatRef] - biasmatRef[fr2>biasmatRef] 
		fr[fr3 < min.delta.frequency] <- NaN
	}

	#select return type
	if(return.type=="ref"){fr[fr<biasmatRef]<- NaN}
	if(return.type=="alt"){fr[fr>biasmatRef]<- NaN}

	#survive stat test max p-value
	if(function.test=="binom.test" & !max.pvalue==1){
		idx <- which(!is.na(fr))
		arn <- arank(x,return.type="names",return.class="matrix") 
		ac <- alleleCounts(x, strand=strand, return.class="array")

		mat1 <- matrix(aperm(ac,c(3,2,1))[aperm(array(matrix(x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,1]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
		],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))

		mat2 <- matrix(aperm(ac,c(3,2,1))[aperm(array(matrix(x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,2]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
		],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))


		allele1 <- mat1[idx] 
		allele2 <- mat2[idx]

		#test the two most expressed alleles
		biasmat1 <- matrix(aperm(mapBias(x, return.class="array"),c(3,2,1))[aperm(array(matrix(
			x@variants, ncol=length(x@variants),
			 nrow=nrow(x), byrow=TRUE)==arn[,1]
		 ,dim=c(nrow(x), length(x@variants), ncol=ncol(x))),c(2,3,1))
		],ncol=ncol(x),nrow=nrow(x), byrow=TRUE, dimnames=list(rownames(x),colnames(x)))


		biasAllele1 <- biasmat1[idx] 
		
		#maybe put a check here later that bias1 and bias2 sum to 1

		ml <- mapply(allele1,allele2,biasAllele1,FUN=function(x,y,z){
			binom.test(c(x,y), p = z)[[3]]
		})

		#cerate matrix in same size as fr	
		pv <- fr	
		pv[idx] <- ml

	}else{stop("function.test must be binom.test")}

	#set non p-value survivors to na
	fr[!(pv<max.pvalue)] <- NaN

	if(return.class=="logical"){
		matrix(!is.na(fr),ncol=ncol(x),nrow=nrow(x), dimnames=list(rownames(x),colnames(x)))
	}else if(return.class=="matrix"){
		stop("return.class matrix as option is not available atm")
	}

})

