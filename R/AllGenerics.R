#Constructor Documentation redirects to sangerseq class documentation 
#' @export 
#' @rdname sangerseq-class
setGeneric("sangerseq", function(obj) standardGeneric("sangerseq"))


#' Make Basecalls
#' 
#' Updates a \code{\link{sangerseq}} class object to contain primary and
#' secondary peak calls.
#' 
#' Scf files do not contain secondary basecalls and ABIF files sometimes (but 
#' not always) contain secondary peak calls that show the secondary peak even if
#' clearly a negative peak. This is problematic in sequence reads where 
#' heterozygous sequence data is contained in the chromatogram data. 
#' \code{makeBaseCalls} identifies basecall windows containing more than one 
#' peak using the provided cutoff ratio and then makes heterozygous calls for 
#' those windows. The primarySeq will always contain the base corresponding to 
#' the maximum peak amplitude within the window. The secondaryPeak will have the
#' same base if the peak was classified as a homozygous peak, the base 
#' corresponding to the second tallest peak if two peaks were above the cutoff, 
#' or an ambiguous base if more than two peaks were identified in the window 
#' that are higher than the cutoff ratio.
#' 
#' @param obj \code{\link{sangerseq}} class object
#' @param ratio cutoff ratio for separating signal and noise. Ratio is relative
#'   to maximum peak in basecall window.
#'   
#' @return \code{\link{sangerseq}} s4 object with new basecalls in the
#' primarySeq and secondarySeq fields. Matrix values are also updated to reflect
#' newly called base positions and amplitudes of all peaks within window.
#' 
#' @seealso \code{\link{chromatogram}}, \code{\link{setAllelePhase}}, 
#' \code{\link{sangerseq}}
#' 
#' @examples
#' hetsangerseq <- readsangerseq(system.file("extdata", 
#'                                           "heterozygous.ab1", 
#'                                           package = "sangerseqR"))
#' hetcalls <- makeBaseCalls(hetsangerseq, ratio = 0.33)
#' 
#' @export

setGeneric("makeBaseCalls", 
           function(obj, ratio=.33) standardGeneric("makeBaseCalls"))

#' Generate Chromatogram
#' 
#' Generates a chromatogram from a \code{\link{sangerseq}} class object.
#' 
#' This function outputs a chromatogram to the current device or to a PDF file 
#' (filename is not NULL). Primary, Secondary or both basecalls can be shown if 
#' they are contained in the \code{\link{sangerseq}} object provided. What is 
#' shown will depend on how they were generated. If generated and provided by 
#' the ABIF or SCF file, then it will show the primary calls on the first line 
#' and the secondary calls on the second line. If generated by 
#' \code{\link{makeBaseCalls}}, then they will show the maximum and alternate 
#' basecalls only for heterozygous peaks. Finally, if the 
#' \code{\link{setAllelePhase}} has been run on the object, then the first line 
#' is the reference sequence and the second line is the alternate allele.
#' 
#' The range of the trace data shown is trimmed to the called sequence by
#' default. Setting \code{trim5} and \code{trim3} to NULL will show the complete
#' trace 5' and 3' of the called bases, respectively. This will generally create
#' a very long trace. Conversely, setting \code{trim5} and \code{trim3} to an
#' integer > 0 will hide the data corresponding to that number of bases at each
#' end. For example, \code{trim5=10} and \code{trim3=20} will remove 10 bases
#' from the 5' end and 20 bases from the 3' end.
#' 
#' Several output parameters can also be set to control how the figure appears.
#' However, it should be noted that if the current device is too small, R will
#' return an error and not show the chromatogram. This is common with long
#' sequence reads. To bypass this error, we recommend that the user set
#' \code{filename}. This will cause the chromatogram to be saved to a PDF file
#' in the current working directory.
#' 
#' @param obj \code{\link{sangerseq}} class object
#' @param trim5 Number of bases to trim from the beginning of the sequence.
#' @param trim3 Number of bases to trim from the end of the sequence.
#' @param showcalls Which basecall sequence to show. Any value other than
#'   "primary", "secondary" or "both" will result in basecalls not being shown.
#' @param width Approximate number of bases per line.
#' @param height Height of each plot row. Ignored by some devices.
#' @param cex.mtext Size factor for the text in the margins.
#' @param cex.base Size factor for the basecall text.
#' @param ylim Peaks larger than this many times the IQR larger than the median
#'   will be cutoff.
#' @param filename Name of PDF file to save to. A "NULL" value outputs the
#'   chromatogram to the current device.
#' @param showtrim If true, highlights trimmed region instead of hiding it.
#' @param showhets Whether or not to highlight heterozygous positions.
#'   
#' @return A plot showing the chromatogram tracedata and, optionally, basecalls.
#' 
#' @seealso \code{\link{makeBaseCalls}}, \code{\link{setAllelePhase}}, 
#' \code{\link{sangerseq}}
#' 
#' @examples
#' hetsangerseq <- readsangerseq(system.file("extdata", 
#'                                           "heterozygous.ab1", 
#'                                           package = "sangerseqR"))
#' hetcalls <- makeBaseCalls(hetsangerseq, ratio = 0.33)
#' chromatogram(hetcalls, width = 100, height = 2, trim5 = 50, trim3 = 100, 
#'              showcalls = "both", filename = "chromatogram.pdf")
#' 
#' @export

setGeneric("chromatogram", 
           function(obj, trim5=0, trim3=0, 
                    showcalls=c("primary", "secondary", "both", "none"), 
                    width=500, height=NA, cex.mtext=1, cex.base=1, ylim=2, 
                    filename=NULL, showtrim=FALSE, showhets=TRUE) 
  standardGeneric("chromatogram")
)

#' Set Reference and Alternate Alleles
#' 
#' Parses the Primary and Secondary Sequences into Reference and Alternate
#' Alleles
#' 
#' When multiple heterozygous basecalls are made, it is generally unclear which
#' calls are in phase with each other. This function takes a reference sequence
#' for one of the alleles to match the primary and secondary basecalls as
#' reference or alternate allele.
#' 
#' @param obj \code{\link{sangerseq}} class object
#' @param refseq DNAString for character string of reference allele sequence.
#' @param trim5 Number of bases to trim from the beginning of the sequence.
#' @param trim3 Number of bases to trim from the end of the sequence.
#'   
#' @return A \code{\link{sangerseq}} object with the Reference Allele in the
#' primarySeq slot and the Alternate Allele in the secondarySeq slot.
#' 
#' @seealso \code{\link{makeBaseCalls}}, \code{\link{chromatogram}}, 
#' \code{\link{sangerseq}}
#' 
#' @examples
#' #Load Sequences
#' hetsangerseq <- readsangerseq(system.file("extdata", 
#'                                           "heterozygous.ab1", 
#'                                           package = "sangerseqR"))
#' homosangerseq <- readsangerseq(system.file("extdata", 
#'                                            "homozygous.scf", 
#'                                            package = "sangerseqR"))
#' #Make calls on heterozygous sequence to be parsed
#' hetcalls <- makeBaseCalls(hetsangerseq, ratio = 0.33)
#' #Need a reference sequence to set phase. Can get from annotation 
#' #(e.g. Refseq) or another sanger sequencing file
#' ref <- subseq(primarySeq(homosangerseq, string = TRUE), 
#'               start = 30, 
#'               width = 500)
#' #Set the phase
#' hetseqalleles <- setAllelePhase(hetcalls, ref, trim5 = 50, trim3 = 100, as.is = FALSE)
#' #Align to compare alleles
#' pa <- pairwiseAlignment(primarySeq(hetseqalleles), 
#'                                    secondarySeq(hetseqalleles), 
#'                                    type = "global-local") 
#' writePairwiseAlignments(pa)
#' 
#' @export

setGeneric("setAllelePhase", 
    function(obj, refSeq, trim5 = 0, trim3 = 0, asIs = FALSE) {
          standardGeneric("setAllelePhase")
    }
)
           

#accessor and setter function generics
#documented under sangerseq class
#' Sangerseq Accessor Functions
#' 
#' Accessor Functions allow the user to retrieve results from and assign values
#' to \code{\link{sangerseq-class}} objects.
#' 
#' @param obj sangerseq object to be manipulated
#' @param string TRUE/FALSE. FALSE (default) returns a \code{\link{DNAString}}
#'   class object. TRUE returns the DNA sequence as a character string.
#' @param value The value to set the slot to.
#'   
#' @seealso \code{\link{sangerseq-class}}
#' 
#' @examples
#' hetsangerseq <- readsangerseq(system.file("extdata", 
#'                                           "heterozygous.ab1", 
#'                                           package = "sangerseqR"))
#' primarySeq(hetsangerseq)
#' secondarySeq(hetsangerseq, string=TRUE)
#' primarySeqID(hetsangerseq)
#' primarySeqID(hetsangerseq) <- "Some string"
#' primarySeqID(hetsangerseq)
#' 
#' @export
#' @docType methods
#' @rdname sangerseqAccessors
#' @aliases primarySeqID

setGeneric("primarySeqID", function(obj) standardGeneric("primarySeqID"))

#' @export
#' @rdname sangerseqAccessors
#' @aliases primarySeqID<-

setGeneric("primarySeqID<-", 
           function(obj, value) standardGeneric("primarySeqID<-"))

#' @export
#' @rdname sangerseqAccessors
#' @aliases primarySeq

setGeneric("primarySeq", 
           function(obj, string=FALSE) standardGeneric("primarySeq"))

#' @export
#' @rdname sangerseqAccessors
#' @aliases primarySeq<-

setGeneric("primarySeq<-", 
           function(obj, value) standardGeneric("primarySeq<-"))

#' @export
#' @rdname sangerseqAccessors
#' @aliases secondarySeqID

setGeneric("secondarySeqID", 
           function(obj) standardGeneric("secondarySeqID"))

#' @export
#' @rdname sangerseqAccessors
#' @aliases secondarySeqID<-

setGeneric("secondarySeqID<-", 
           function(obj, value) standardGeneric("secondarySeqID<-"))

#' @export
#' @rdname sangerseqAccessors
#' @aliases secondarySeq

setGeneric("secondarySeq", 
           function(obj, string=FALSE) standardGeneric("secondarySeq"))

#' @export
#' @rdname sangerseqAccessors
#' @aliases secondarySeq<-

setGeneric("secondarySeq<-", 
           function(obj, value) standardGeneric("secondarySeq<-"))

#' @export
#' @rdname sangerseqAccessors
#' @aliases traceMatrix

setGeneric("traceMatrix", 
           function(obj) standardGeneric("traceMatrix"))

#' @export
#' @rdname sangerseqAccessors
#' @aliases traceMatrix<-

setGeneric("traceMatrix<-", 
           function(obj, value) standardGeneric("traceMatrix<-"))

#' @export
#' @rdname sangerseqAccessors
#' @aliases peakPosMatrix

setGeneric("peakPosMatrix", 
           function(obj) standardGeneric("peakPosMatrix"))

#' @export
#' @rdname sangerseqAccessors
#' @aliases peakPosMatrix<-

setGeneric("peakPosMatrix<-", 
           function(obj, value) standardGeneric("peakPosMatrix<-"))

#' @export
#' @rdname sangerseqAccessors
#' @aliases peakAmpMatrix

setGeneric("peakAmpMatrix", 
           function(obj) standardGeneric("peakAmpMatrix"))

#' @export
#' @rdname sangerseqAccessors
#' @aliases peakAmpMatrix<-

setGeneric("peakAmpMatrix<-", 
           function(obj, value) standardGeneric("peakAmpMatrix<-"))

