\name{alignSequences}
\alias{alignSequences}
\title{
Align sequences via clustalo
}
\description{
Align sequences for an object of class LowMACA
}
\usage{
alignSequences(object, clustalo_filename=NULL , mail=NULL , 
	perlCommand="perl", use_hmm=FALSE, datum=FALSE)
}
\arguments{
  \item{object}{ an object of class LowMACA containing at least 2 sequences.}
  \item{clustalo_filename}{ a character string that contains the file name where 
  			clustal omega alignment file will be stored. In case it's NULL
  			no file will be written. Default=NULL}
  \item{mail}{ a character string indicating the email address 
  				where error report should be sent in web mode}
  \item{perlCommand}{ a character string containing the path to Perl executable. 
				if missing, "perl" will be used as default}
  \item{use_hmm}{When analysing Pfam sequences, it is possible to
				use the Hidden Markov Model (HMM) of the 
				specific Pfam to align the sequences.
				Default is FALSE.}
  \item{datum}{When analysing Pfam sequences, use all the genes
				that belong to the Pfam to generate the
				alignment. This creates a unique mapping between
				individual residues and consensus sequence, 
				disregarding the set of sequences that are 
				selected for the analysis.
				Default is FALSE.}
}
\details{
This method launches a system call to clustalo aligner
and optionally creates a fasta file in clustal format.
A warning is returned if at least one sequence has a pairwise similarity below 
20\% to any other sequence.
If only one sequence is passed to \code{alignSequences}, the alignment will be skipped, 
but no warning will be raised.
If mail is not NULL, a local installation of clustal omega 
is no longer required and the alignment is 
performed using clustal omega EBI web service. A limit of 2000 sequences is set 
in this case and perl must be installed in the system
}
\value{
The method returns an object of class LowMACA updating the slot \code{alignment}.
See \code{\link{lmAlignment}}
}
\references{
\href{http://www.ncbi.nlm.nih.gov/pubmed/12112692}{Trident Score}
\href{http://www.clustal.org/omega/}{Clustal Omega}
\href{http://www.ebi.ac.uk/Tools/webservices/services/msa/clustalo_soap}{Clustal Omega Web Service}
}
\author{
Stefano de Pretis, Giorgio Melloni
}
\section{Warning}{
When a sequence has a similarity below 20\%, a warning is raised. In order to produce strong results 
in terms of conservation of multiple mutations, consider to remove that sequence 
from the analysis. The alignment will obviously change.
}
\seealso{
\code{\link{getMutations}} , \code{\link{mapMutations}} , \code{\link{setup}}
}
\examples{
#Create an object of class LowMACA for RAS domain family
lm <- newLowMACA(pfam="PF00071" , genes=c("KRAS" , "NRAS" , "HRAS"))
#Align sequences using local installation of clustalo
lm <- alignSequences(lm)
#Web service clustalomega aligner
lm <- alignSequences(lm , mail="lowmaca@gmail.com")
#Use HMM to align
lm <- alignSequences(lm , use_hmm=TRUE)
#Use "datum"
lm <- alignSequences(lm , datum=TRUE)
}
