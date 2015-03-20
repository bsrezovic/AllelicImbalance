#' A package meant to provide all basic functions for high-throughput allele
#' specific expression analysis
#' 
#' Package \code{AllelicImbalance} has functions for importing, filtering and
#' plotting high-throughput data to make an allele specific expression
#' analysis. A major aim of this package is to provide functions to collect as
#' much information as possible from regions of choice, and to be able to
#' explore the allelic expression of that region in detail.
#' 
#' \tabular{ll}{ Package: \tab AllelicImbalance\cr Type: \tab Package\cr
#' Version: \tab 1.2.0\cr Date: \tab 2014-08-24\cr License: \tab GPL-3\cr }
#' 
#' @name AllelicImbalance-package
#' @aliases AllelicImbalance-package AllelicImbalance
#' @docType package
#' @section Overview - standard procedure: Start out creating a \code{GRange}
#' object defining the region of interest. This can also be done using
#' \code{getAreaFromGeneNames} providing gene names as arguments. Then use
#' \code{BamImpGAList} to import reads from that reagion and find potential
#' SNPs using \code{scanForHeterozygotes}. Then retrieve the allele counts of
#' heterozygote sites by the function \code{getAlleleCount}. With this data
#' create an \code{ASEset}. At this point all pre-requisites for a 'basic'
#' allele specific expression analysis is available. Two ways to go on could be
#' to apply \code{\link{chisq.test}} or \code{\link{barplot}} on this ASEset
#' object.
#' @author Author: Jesper Robert Gadin Author: Lasse Folkersen
#' 
#' Maintainer: Jesper Robert Gadin <j.r.gadin@@gmail.com>
#' @seealso \itemize{ \item code?ASEset }
#' @references Reference to published application note (work in progress)
#' @keywords package
#'
#' @import methods
#' @import BiocGenerics
#' @import S4Vectors
#' @import grid
#' @import GenomeInfoDb
#' 
#' 
#' @importFrom IRanges IRanges
#' @importFrom IRanges ranges
#' @importFrom IRanges subsetByOverlaps
#' @importFrom IRanges coverage
#' @importFrom IRanges unlist
#' @importFrom IRanges which
#' @importFrom IRanges resize
#' @importFrom IRanges IntervalTree
#' 
# @importFrom S4Vectors DataFrame
# @importFrom S4Vectors SimpleList
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges granges
#' @importFrom GenomicRanges GRangesList
#' @importFrom GenomicRanges unlist
#' @importFrom GenomicRanges reduce
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges rowData
#' @importFrom GenomicRanges rowRanges
#' @importFrom GenomicRanges colData
#' @importFrom GenomicRanges SummarizedExperiment
#' @importFrom GenomicRanges grglist
#' @importFrom GenomicRanges flank
#' @importFrom GenomicRanges assays
#' @importFrom GenomicRanges assays<-
#' 
#' @importFrom Biostrings subseq
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings append
#' 
#' @importFrom Rsamtools scanBamWhat
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom Rsamtools scanBamHeader
#' @importFrom Rsamtools indexBam
#' @importFrom Rsamtools BamFileList
#' @importFrom Rsamtools scanBam
#' @importFrom Rsamtools ScanBcfParam
#' @importFrom Rsamtools scanBcf
#' @importFrom Rsamtools BcfFile
#' @importFrom Rsamtools BcfFileList
#' @importFrom Rsamtools indexBcf
#' @importFrom Rsamtools indexFa
#' @importFrom Rsamtools FaFile
#' @importFrom Rsamtools scanFa
#' @importFrom Rsamtools scanFaIndex
#' @importFrom Rsamtools applyPileups
#' @importFrom Rsamtools ApplyPileupsParam
#' @importFrom Rsamtools PileupFiles
#' 
#' @importClassesFrom GenomicAlignments GAlignments
#' @importClassesFrom GenomicAlignments GAlignmentPairs
#' @importClassesFrom GenomicAlignments GAlignmentsList
#' 
#' @importFrom GenomicAlignments pileLettersAt
#' @importFrom GenomicAlignments cigar
#' @importFrom GenomicAlignments cigarToRleList
#' @importFrom GenomicAlignments GAlignments
#' @importFrom GenomicAlignments GAlignmentPairs
#' @importFrom GenomicAlignments GAlignmentsList
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicAlignments readGAlignmentPairs
#' @importFrom GenomicAlignments cigarWidthAlongReferenceSpace
#' 
# @importFrom VariantAnnotation readVcf
# @export readVcf
# @importFrom VariantAnnotation readGT
# @export readGT
#'
#' @importFrom lattice barchart
#' @importFrom lattice xyplot
#' @importFrom lattice trellis.par.set
#' @importFrom lattice panel.abline
#' @importFrom lattice panel.text
#' @importFrom lattice panel.barchart
# @importFrom gridExtra grid.arrange
#'
#' @importFrom seqinr write.fasta
#'
#' @importClassesFrom GenomicRanges GRangesList
#' @importClassesFrom GenomicRanges GRanges
# @importClassesFrom Gviz DataTrack
# @importClassesFrom Gviz AnnotationTrack
# @importClassesFrom Gviz DetailsAnnotationTrack
# @importClassesFrom Gviz GeneRegionTrack
# @importClassesFrom Gviz GenomeAxisTrack
#
# @importFrom Gviz DataTrack
# @importFrom Gviz AnnotationTrack
# @importFrom Gviz DetailsAnnotationTrack
# @importFrom Gviz GeneRegionTrack
# @importFrom Gviz GenomeAxisTrack
#
# @importFrom Gviz plotTracks
# @export plotTracks
#' 
#' @importClassesFrom GenomicFeatures TxDb
#' @importClassesFrom AnnotationDbi AnnotationDb
#' @importClassesFrom AnnotationDbi OrgDb
#' @importClassesFrom AnnotationDbi ChipDb
#' @importClassesFrom AnnotationDbi InparanoidDb
#' @importClassesFrom AnnotationDbi GODb
#' @importClassesFrom AnnotationDbi ReactomeDb

#' @importMethodsFrom AnnotationDbi columns
#' 
#' @importFrom AnnotationDbi select
#' @importFrom AnnotationDbi keys
#' 
#' @importFrom GenomicFeatures exons
#' @importFrom GenomicFeatures transcripts
#' @importFrom GenomicFeatures cds
#' @importFrom GenomicFeatures isActiveSeq
#' @importFrom GenomicFeatures 'isActiveSeq<-'
NULL

## @importFrom SNPlocs.Hsapiens.dbSNP.20120608 getSNPlocs

#' GRvariants object
#' 
#' this data is a \code{GRanges} object that contains the ranges for three
#' example SNPs.
#' 
#' 
#' @name GRvariants
#' @docType data
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{reads}} which is another example
#' object }
#' @keywords data object example
#' @examples
#' 
#' #load example data
#' data(GRvariants)
#' 
#' 
NULL


#' reads object
#' 
#' This data set corresponds to the BAM-file data import illustrated in the
#' vignette. The data set consists of a chromosome 17 region from 20 RNA-seq
#' experiments of HapMap samples.
#' 
#' 
#' @name reads
#' @docType data
#' @author Jesper R. Gadin, Lasse Folkersen
#' @seealso \itemize{ \item The \code{\link{GRvariants}} which is another
#' example object }
#' @references Montgomery SB et al. Transcriptome genetics using second
#' generation sequencing in a Caucasian population. Nature. 2010 Apr
#' 1;464(7289):773-7.
#' @keywords data object example
#' @examples
#' 
#' ##load eample data (Not Run)  
#' #data(reads)
#' 
NULL 


#' ASEset.sim object
#' 
#' ASEset with simulated data with SNPs within the first 200bp of chromosome 17,
#' which is required to have example data for the refAllele function.
#' 
#' @name ASEset.sim
#' @docType data
#' @author Jesper R. Gadin, Lasse Folkersen
#' @keywords data object example
#' @examples
#' 
#' ##load eample data (Not Run)  
#' #data(ASEset.sim)
#' 
NULL 
