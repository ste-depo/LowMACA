\name{LowMACA-class}
\Rdversion{0.99.5}
\docType{class}
\alias{LowMACA-class}
\alias{alignSequences,LowMACA-method}
\alias{bpAll,LowMACA-method}
\alias{entropy,LowMACA-method}
\alias{getMutations,LowMACA-method}
\alias{lfm,LowMACA-method}
\alias{lmPlot,LowMACA-method}
\alias{mapMutations,LowMACA-method}
\alias{nullProfile,LowMACA-method}
\alias{parallelize,LowMACA-method}
\alias{parallelize<-,LowMACA-method}
\alias{lmAlignment,LowMACA-method}
\alias{lmMutations,LowMACA-method}
\alias{lmEntropy,LowMACA-method}
\alias{lmParams,LowMACA-method}
\alias{lmParams<-,LowMACA-method}
\alias{protter,LowMACA-method}
\alias{setup,LowMACA-method}
\alias{show,LowMACA-method}
\alias{lfmSingleSequence,LowMACA-method}
\alias{lmPlotSingleSequence,LowMACA-method}
\title{Class \code{"LowMACA"}}
\description{
LowMACA class object describing the properties of mutations mapped on pfam domains or proteins
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{newLowMACA(genes, pfam)}.
}
\section{Constructor}{
newLowMACA(genes=character_vector , pfam=character_vector)
}
\section{Slots}{
\describe{
  \item{arguments}{
    Object of class \code{"list"} with 6 elements:
    \itemize{
        \item{genes}{ : vector of selected genes for the analysis in Hugo names format. NULL if \code{mode="pfam"}.}
        \item{pfam}{ : vector of selected domains for the analysis in pfam ids format. NULL if \code{mode="genes"}.}
        \item{input}{ : data.frame describing the input data as gene symbols, pfam ids, entrez ids, 
                            envelope start and end of the domain relative to the 
                            protein, name of the canonical protein in uniprot format, amino acidic sequence.}
        \item{mode}{ : character. automatically set by the constructor as either \code{"pfam"} or \code{"genes"}.
                            If \code{pfam=NULL} then \code{mode="genes"}, \code{"pfam"} otherwise.}
        \item{params}{ : named list of starting parameters for the LowMaca analysis. Call \code{lmParams(object)} to show default. 
                            See \code{\link{lmParams}} for further details.}
        \item{parallelize}{ : named list of logicals. \code{getMutations}=FALSE is the default for the \code{\link{getMutations}} 
                                    method and \code{makeAlignment}=TRUE is the default for the \code{\link{alignSequences}} method.
                                    See \code{\link{parallelize}} for further details.}
        }
    }
  \item{alignment}{
    Object of class \code{"list"} with 4 elements:
      \itemize{
        \item{ALIGNMENT}{ : data.frame of the result of the alignment. 
                        Every row represents a position of a sequence and the relative mapping to the consensus sequence.}
        \item{SCORE}{ : list of two elements. \code{DIST_MAT} is a matrix of pairwise similarities
                              between sequences as described by clustalo. \code{SUMMARY_SCORE} is a dataframe 
                              of summary descriptive statistics of the \code{DIST_MAT} matrix}
        \item{CLUSTAL}{ : an object of class \code{\link[Biostrings]{AAMultipleAlignment}} from package Biostrings}
        \item{df}{ : a data.frame describing the consensus sequence, its per-position degree 
                          of conservation and its mutations null profile density.
                          See \code{\link{entropy}} and \code{\link{lmPlot}} for further details}
        }
  }
  \item{mutations}{
      Object of class \code{"list"} with 3 elements:
      \itemize{
        \item{data}{ : data.frame derived from the query to the cBioPortal query, \code{\link[cgdsr]{getMutationData}}
                            Every row represents a mutation stratified by position, gene and tumor type.}
        \item{freq}{ : data.frame of absolute frequency of mutation stratified by gene and tumor type.}
        \item{aligned}{ : matrix representing the number of mutations at every position in the consensus sequence 
                                (columns) and in each original sequence (rows)}
        }
    }
  \item{entropy}{
    Object of class \code{"list"} with 5 elements:
      \itemize{
        \item{bw}{ : numeric value. user defined bandwidth for the function \code{\link{entropy}}}
        \item{uniform}{ : function that generate the uniform null profile}
        \item{absval}{ : numeric value. Shannon entropy of the mutation data profile according to the defined bandwidth}
        \item{log10pval}{ : numeric value. pvalue of the entropy test in -log10 scale}
        \item{pvalue}{ : numeric value. pvalue of the entropy test}
      }
    }
  }
}
\section{Methods}{
  \describe{
    \item{alignSequences}{\code{alignSequences(object = "LowMACA")}: ... }
    \item{bpAll}{\code{bpAll(object = "LowMACA")}: ... }
    \item{entropy}{\code{entropy(object = "LowMACA")}: ... }
    \item{getMutations}{\code{getMutations(object = "LowMACA")}: ... }
    \item{lfm}{\code{lfm(object = "LowMACA")}: ... }
    \item{lmPlot}{\code{lmPlot(object = "LowMACA")}: ... }
    \item{mapMutations}{\code{mapMutations(object = "LowMACA")}: ... }
    \item{nullProfile}{\code{signature(object = "LowMACA")}: ... }
    \item{parallelize}{\code{parallelize(object = "LowMACA")}: ... }
    \item{parallelize<-}{\code{signature(object = "LowMACA")}: ... }
    \item{lmParams}{\code{params(x = "LowMACA")}: ... }
    \item{lmParams<-}{\code{signature(object = "LowMACA")}: ... }
    \item{protter}{\code{protter(object = "LowMACA")}: ... }
    \item{setup}{\code{setup(object = "LowMACA")}: ... }
    \item{show}{\code{show(object = "LowMACA")}: ... }
    \item{lfmSingleSequence}{\code{lfmSingleSequence(object = "LowMACA")}: ... }
    \item{lmPlotSingleSequence}{\code{lmPlotSingleSequence(object = "LowMACA")}: ... }
	 }
}
\references{
\href{https://cgsb.genomics.iit.it/wiki/projects/LowMACA}{LowMACA website}
}
\author{
Stefano de Pretis, Giorgio Melloni
}
\seealso{
\code{\link{newLowMACA}}
}
\examples{
#ANALYSIS OF SOME OF THE PROTEINS THAT SHARE THE HOMEOBOX DOMAIN
#Genes to analyze
Genes <- c("ADNP","ALX1","ALX4","ARGFX","CDX4","CRX"
  		,"CUX1","CUX2","DBX2","DLX5","DMBX1","DRGX"
			,"DUXA","ESX1","EVX2","HDX","HLX","HNF1A"
			,"HOXA1","HOXA2","HOXA3","HOXA5","HOXB1","HOXB3"
			,"HOXD3","ISL1","ISX","LHX8")
#Pfam to analyze
Pfam <- "PF00046"
#Construct a new LowMACA object
lm <- newLowMACA(genes=Genes , pfam=Pfam)
#Change some parameters
lmParams(lm)[['tumor_type']] <- c("skcm" , "stad" , "ucec" , "luad" , "lusc" , "coadread" , "brca")
lmParams(lm)[['min_mutation_number']] <- 1
lmParams(lm)[['density_bw']] <- 0
#Run if you have clustalo installed
lm <- setup(lm)
#Calculate staistics
lm <- entropy(lm)
#Retrieve original mutations
lfm(lm)
#Plot
bpAll(lm)
lmPlot(lm)
protter(lm)
}
\keyword{classes}
