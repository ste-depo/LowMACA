\name{lmParams}
\alias{lmParams}
\alias{lmParams<-}
\title{
Show and set parameters
}
\description{
Method for objects of class LowMACA. It can show the most important user-definable parameters
for a LowMACA analysis and allows to change them.
}
\usage{
lmParams(object)
lmParams(object) <- value
}
\arguments{
\item{object}{ an object of class LowMaca}
\item{value}{ a named list containing:
  \enumerate{
    \item \code{mutation_type} a character string among: 'missense' , 'truncating' , 'silent' ,'all'. Default 'missense'
    \item \code{tumor_type} a character vector or string containing the tumor type barcode of the data in cBioPortal.
                              Default 'all'.
    \item \code{min_mutation_number} an integer value describing the minimum number of mutations 
                                    accepted for a sequence. If a sequence does not harbor a sufficient 
                                    number of mutations is discarded from the analysis. Default is 1
    \item \code{density_bw} either a numeric value or 'auto'. A numeric value is passed directly 
                              to the function density while putting 0 will not launch density at all (every position 
                              is not aggregated to the surrounded ones). 'auto' will let the simulation 
                              decide according to the Silverman's rule of thumb the correct bandwidth.
                              Default is 0.
    \item \code{clustal_cmd} path to clustalo executable
    \item \code{use_hmm} When analysing Pfam sequences, it is possible to
                              use the Hidden Markov Model (HMM) of the 
                              specific Pfam to align the sequences.
                              Default is FALSE.
    \item \code{datum} When analysing Pfam sequences, use all the genes
                              that belong to the Pfam to generate the
                              alignment. This creates a unique mapping between
                              individual residues and consensus sequence, 
                              disregarding the set of sequences that are 
                              selected for the analysis.
                              Default is FALSE.
      }
    }
}
\details{
LowMACA is a suite of tool that analyze conserved mutations, 
so it looks for clusters of gain of function alterations.
With 'missense' mutation_type we intend all those mutations that change 
the original DNA but do not create stop codon nor alter the reading frame 
(these mutations are collectively defined as 'truncating' mutations).
In addition we let the possibility to also choose 'silent' mutations 
even though they are currently not supported by the cBioPortal.
To see all the available tumor types to run a LowMACA analysis, simply run \code{showTumorType}.
The parameter \code{density_bw} has a strong effect on the statistical analysis of LowMACA. With the default bandwidth (0), 
the Shannon entropy calculation becomes descrete, while the continuos version is used in all the other cases.
}
\value{
If \code{lmParams} is used as a show method it returns a named list of 5 elements:
\code{mutation_type='missense'} , \code{tumor_type='all'} , 
\code{min_mutation_number=1} , \code{density_bw=0} , \code{clustal_cmd='clustalo'}
}
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{showTumorType}} \code{\link{getMutations}} \code{\link{entropy}} \code{\link[stats]{density}}
}
\examples{
#Construct a LowMACA object
lm <- newLowMACA(pfam="PF12906")
#Show default parameters
lmParams(lm)
#Change all parameters
lmParams(lm) <- list(mutation_type='all' 
                    , tumor_type=c('skcm','brca') 
                    , min_mutation_number=0 
                    , density_bw=0
                    , clustal_cmd='clustalo'
                    , use_hmm=FALSE
                    , datum=FALSE)
#Change just one parameter
lmParams(lm)[['tumor_type']] <- 'prad'
}
