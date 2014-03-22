\name{funGetAnnotation-functions}
%annotationDb
\alias{getGenesFromAnnotation}
\alias{getGenesVector}
%transcriptDb
\alias{getExonsFromAnnotation}
\alias{getTranscriptsFromAnnotation}
\alias{getCDSFromAnnotation}

\alias{getExonsVector}
\alias{getTranscriptsVector}
\alias{getCDSVector}

%common wrapper
\alias{getAnnotationDataFrame}

\title{AnnotationDb wrappers}
\description{
    These functions acts as wrappers to retrieve information from annotation database objects (\code{annotationDb objects}) or (\code{transcriptDb objects})
}

\usage{
	getGenesFromAnnotation(OrgDb,GR,leftFlank=0,rightFlank=0,
    getUCSC=FALSE,verbose=FALSE)
  getExonsFromAnnotation(TxDb,GR,leftFlank=0,rightFlank=0,
      verbose=FALSE)
  getTranscriptsFromAnnotation(TxDb,GR,leftFlank=0,
      rightFlank=0,verbose=FALSE)
  getCDSFromAnnotation(TxDb,GR,leftFlank=0,rightFlank=0,
    verbose=FALSE)

  getGenesVector(OrgDb, GR, leftFlank=0, rightFlank=0, verbose=FALSE)
  getExonsVector(TxDb, GR,leftFlank=0,rightFlank=0, verbose=FALSE)
  getTranscriptsVector(TxDb, GR,leftFlank=0,rightFlank=0, verbose=FALSE)
  getCDSVector(TxDb, GR,leftFlank=0,rightFlank=0,verbose=FALSE)

  getAnnotationDataFrame(GR,strand="+",annotationType=NULL,
    OrgDb=NULL, TxDb=NULL,verbose=FALSE)

}

\arguments{
  \item{OrgDb}{An \code{OrgDb} object}
  \item{GR}{A \code{GenomicRanges} object with sample area}
  \item{leftFlank}{An \code{integer} specifying number of additional nucleotides around the SNPs for the leftFlank}
  \item{rightFlank}{An \code{integer} specifying number of additional nucleotides around the SNPs for the rightFlank}
  \item{getUCSC}{A \code{logical} indicating if UCSC transcript IDs should also be retrieved}
  \item{TxDb}{A \code{transcriptDb} object}
  \item{strand}{Two options,"+" or "-"}
  \item{annotationType}{select one or more from "gene", "exon", "transcript", "cds". }
  \item{verbose}{A \code{logical} making the functions more talkative}

}
\details{
	These functions retrieve regional annotation from OrgDb or TxDb objects, when given GRanges objects. 
}
\value{
    %OrgDb
  The \code{getGenesFromAnnotation}{
    	 function will return a \code{GRanges object} with ranges over the genes in the region.}
  
   The \code{getGenesVector}{
       function will return a character vector where each element are gene symbols separated by comma}

    %transcriptDb
   The \code{getExonsFromAnnotation}{
       function will return a \code{GRanges object} with ranges over the exons in the region.}
  
   The \code{getTranscriptsFromAnnotation}{
       function will return a \code{GRanges object} with ranges over the transcripts in the region.} 
  
   The \code{getCDSFromAnnotation}{
       function will return a \code{GRanges object} with ranges over the CDSFs in the region.}  

   The \code{getExonsVector}{
       function will return a character vector where each element are exons separated by comma}
  
   The \code{getTranscriptsVector}{
       function will return a character vector where each element are transcripts separated by comma}
  
   The \code{getCDSVector}{
       function will return a character vector where each element are CDSs separated by comma}

    The \code{getAnnotationDataFrame}{
       function will return a data.frame with annotations. This function is used internally by i.e. the barplot-function}

}
\author{Jesper R. Gadin, Lasse Folkersen}

\examples{

  data(ASEset)
  require(org.Hs.eg.db)
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  OrgDb <- org.Hs.eg.db
  TxDb <- TxDb.Hsapiens.UCSC.hg19.knownGene

  #use for example BcfFiles as the source for SNPs of interest
  GR <- rowData(ASEset)
  #get annotation
  g <- getGenesFromAnnotation(OrgDb,GR)
  e <- getExonsFromAnnotation(TxDb,GR)
  t <- getTranscriptsFromAnnotation(TxDb,GR)
  c <- getCDSFromAnnotation(TxDb,GR)
  
}

\keyword{ genes }
\keyword{ exons }
\keyword{ transcripts }
\keyword{ CDS }
\keyword{ annotation }


