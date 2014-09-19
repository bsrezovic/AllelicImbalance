
#' inference of SNPs of ASEset
#' 
#' inference of SNPs
#' 
#' threshold.frequency is the least fraction needed to classify as bi tri or
#' quad allelic SNPs. If 'all' then all of bi tri and quad allelic SNPs will use the same
#' threshold. Everything under the treshold will be regarded as noise. 'all' will return 
#' a matrix with snps as rows and uni bi tri and quad will be columns. For this function
#' Anything that will return TRUE for tri-allelicwill also return TRUE for uni and bi-allelic
#' for the same SNP an Sample.
#' 
#'
#' @name inferAlleles
#' @rdname inferAlleles
#' @aliases inferAlleles
#' inferAlleles,ASEset-method 
#' @docType methods
#' @param x ASEset
#' @param strand strand to infer from
#' @param inferOver 'eachSample' or 'allSamples' 
#' @param allow.NA treat NA as zero when TRUE
#' @param return.type 'uni' 'bi' 'tri' 'quad' 'all' 
#' @param threshold.count.sample least amount of counts to try to infer allele
#' @param threshold.frequency least fraction to classify (see details)
#' @author Jesper R. Gadin
#' @keywords infer
#' @examples
#' 
#' data(ASEset)
#' i <- inferAlleles(ASEset)
#' 
#' @exportMethod inferAlleles

setGeneric("inferAlleles", function(x,strand="*",return.type="all",
	threshold.frequency=0.2, threshold.count.sample=5,
	inferOver="eachSample", allow.NA=FALSE
	){standardGeneric("inferAlleles")
})
setMethod("inferAlleles", signature(x = "ASEset"), function(x,strand="*",return.type="all",
	threshold.frequency=0.0, threshold.count.sample=0,
	inferOver="eachSample", allow.NA=FALSE
	){

	#find SNP types
	fr <- frequency(x,strand=strand, 
		threshold.count.sample=threshold.count.sample,
		return.class="array")
	
	#apply over alleles
	alleles <- apply(fr,c(1,2), function(x){
		sum(x >= threshold.frequency)
	})

	#disregard samples with NA
	if(!allow.NA){
		alleles[is.na(alleles)] <- 0
	}

	if(inferOver=="eachSample"){
		
		if(return.type%in% c("all","uni")){
			tfuni <- alleles == 1 
		}
		if(return.type%in% c("all","bi")){
			tfbi <- alleles == 2 
		}
		if(return.type%in% c("all","tri")){
			tftri <- alleles == 3 
		}
		if(return.type%in% c("all","quad")){
			tfquad <- alleles == 4 
		}
		
		if(return.type=="all"){
			tfall <- array(c(
				 tfuni | tfbi | tftri | tfquad,
				 tfbi | tftri | tfquad,
				 tftri | tfquad,
				 tfquad),
				 dim=c(nrow(x),ncol(x),4),
				 dimnames=list(rownames(x), colnames(x), c("uni","bi","tri","quad")))
			
		}	

	}
	if(inferOver=="allSamples"){
		#apply over samples
		if(return.type%in% c("all","uni")){
			tfuni <- apply(alleles,1, function(x){
				sum(x == 1 ) >= 1
			})
		}
		if(return.type%in% c("all","bi")){
			tfbi <- apply(alleles,1, function(x){
				sum(x == 2 ) >= 1
			})
		}
		if(return.type%in% c("all","tri")){
			tftri <- apply(alleles,1, function(x){
				sum(x == 3 ) >= 1
			})
		}
		if(return.type%in% c("all","quad")){
			tfquad <- apply(alleles,1, function(x){
				sum(x == 4 ) >= 1
			})
		}
		
		if(return.type=="uni"){
			tfuni
		}else if(return.type=="bi"){
			tfbi
		}else if(return.type=="tri"){
			tftri
		}else if(return.type=="quad"){
			tfquad
		}else if(return.type=="all"){
			tfall <- t(matrix(c(
				 tfuni | tfbi | tftri | tfquad,
				 tfbi | tftri | tfquad,
				 tftri | tfquad,
				 tfquad
				 ), nrow=4,byrow=TRUE,
				 dimnames=list(c("uni","bi","tri","quad"),rownames(x))))
		}	
	}

	if(return.type=="uni"){
		tfuni
	}else if(return.type=="bi"){
		tfbi
	}else if(return.type=="tri"){
		tftri
	}else if(return.type=="quad"){
		tfquad
	}else if(return.type=="all"){
		tfall 
			
	}	
	
})


#' infererence of genotypes from ASEset count data
#' 
#' inference of genotypes
#' 
#' Oftern necessary information to link AI to SNPs outside coding region
#' 
#' @name inferGenotypes
#' @rdname inferGenotypes
#' @aliases inferGenotypes inferGenotypes,ASEset-method
#' @param x ASEset
#' @param strand strand to infer from
#' @param threshold.count.sample least amount of counts to try to infer allele
#' @param threshold.frequency least fraction to classify (see details)
#' @param return.allele.allowed vector with 'bi' 'tri' or 'quad'.
#' 'uni' Always gets returned
#' @param return.class 'matrix' or 'vector'
#' @author Jesper R. Gadin
#' @keywords infer
#' @examples
#' 
#' data(ASEset)
#' g <- inferGenotypes(ASEset)
#' 
#' @exportMethod inferGenotypes

setGeneric("inferGenotypes", function(x, strand="*", return.class="matrix",
	return.allele.allowed = c("bi","tri","quad"),
	threshold.frequency = 0.05, threshold.count.sample = 0
	){ standardGeneric("inferGenotypes")})
					   
setMethod("inferGenotypes", signature(x = "ASEset"), function(x, strand="*", return.class="matrix",
	return.allele.allowed = c("bi","tri","quad"),
	threshold.frequency = 0.05, threshold.count.sample = 0
	){ 
	rn <- arank(x, strand=strand,return.type="names", return.class="matrix")			

	if(return.class=="matrix"){
		al <- inferAlleles(x, strand=strand, threshold.frequency=threshold.frequency,
			threshold.count.sample=threshold.count.sample, return.type="all")	

		ar.rn <- aperm(array(rn, dim=c(nrow(x),4, ncol(x)),
			dimnames=list(rownames(x), c("uni","bi","tri","quad"), colnames(x))), c(1,3,2))
			
		ar.rn[!al] <- ""

		#remove all alleletypes not to be returned
		if(!("quad" %in% return.allele.allowed)){
			ar.rn[,,4] <- ""
		}
		if(!("tri" %in% return.allele.allowed)){
			ar.rn[,,3] <- ""
		}
		if(!("bi" %in% return.allele.allowed)){
			ar.rn[,,2] <- ""
		}

		ret <- apply(ar.rn,c(1,2),function(x){
			tf <- x ==""

			if(sum(!tf)==1){
				paste(x[1],x[1],sep="/")
			}else if(sum(!tf)>1){
				paste(x[!tf],collapse="/")
			}else{
				NA
			}
		})

	} else if(return.class=="vector"){
		#infer allele types
		al <- inferAlleles(x, strand=strand, threshold.frequency=threshold.frequency,
			threshold.count.sample=threshold.count.sample, inferOver="allSamples",
			return.type="all")	

		#remove all not
		rn[!al] <- ""

		#remove all alleletypes not to be returned
		if(!("quad" %in% return.allele.allowed)){
			rn[,4] <- ""
		}
		if(!("tri" %in% return.allele.allowed)){
			rn[,3] <- ""
		}
		if(!("bi" %in% return.allele.allowed)){
			rn[,2] <- ""
		}

		ret <- apply(rn,1,function(x){
			tf <- x ==""
			if(sum(!tf)==1){
				paste(x[1],x[1],sep="/")
			}else{
				paste(x[!tf],collapse="/")
			}
		})
	}
	ret
})

