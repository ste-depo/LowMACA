\name{setup}
\alias{setup}
\title{
Setup of a LowMACA object
}
\description{
A wrapper around \code{alignSequences} , \code{getMutations} 
and \code{mapMutations} in order to execute all these methods at once.
}
\usage{
setup(object, repos = NULL, clustalo_filename=NULL 
	, mail=NULL , perlCommand="perl", use_hmm=FALSE, datum=FALSE)
}
\arguments{
\item{object}{ an object of class LowMACA}
\item{repos}{ a data.frame containing mutations for the specified genes 
            in the LowMACA object in case of custom mutation data. Default NULL}
\item{clustalo_filename}{ a character string that contains the file name where 
			clustal omega alignment file will be stored. In case it's NULL
			no file will be written. Default=NULL}
\item{mail}{ a character string indicating the email address 
  				where error report should be sent in web mode. Default is NULL, 
  				to use a local clustalo installation}
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
If mail is not NULL, a local installation of clustal omega is no longer 
required and the alignment is performed using clustal omega EBI web service. 
A limit of 2000 sequences is set in this case and perl is required with 
XML::Simple and LWP modules installed
}
\value{
An object of class LowMACA with all the updates provided 
by \code{alignSequences} , \code{getMutations} 
and \code{mapMutations} methods.
}
\references{
\href{http://www.ncbi.nlm.nih.gov/pubmed/12112692}{Trident Score}
\href{http://www.clustal.org/omega/}{Clustal Omega}
\href{http://www.ebi.ac.uk/Tools/webservices/services/msa/clustalo_soap}{Clustal Omega Web Service}
}
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{alignSequences}} \code{\link{getMutations}} \code{\link{mapMutations}}
}
\examples{
#Create an object of class LowMACA for RAS domain family
lm <- newLowMACA(pfam="PF00071" , genes=c("KRAS" , "NRAS" , "HRAS"))
#Select a few tumor types
lmParams(lm)$tumor_type <- c("skcm" , "brca" , "coadread")
#Align sequences, get mutation data and map them on consensus
lm <- setup(lm)
#Same as above, but using web service
lm <- setup(lm , mail="lowmaca@gmail.com")
#Use HMM to align
lm <- setup(lm , use_hmm=TRUE)
#Use "datum"
lm <- setup(lm , datum=TRUE)
}
