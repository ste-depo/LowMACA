\name{lfmSingleSequence}
\alias{lfmSingleSequence}
\title{
Show significant clusters of mutations of every 
gene in a LowMACA object without alignment
}
\description{
The method \code{lfmSingleSequence} 
(low frequency mutations in Single Sequence) launch \code{lfm} method
on every gene or domain inside a LowMACA object without aligning the sequences
}
\usage{
lfmSingleSequence(object , metric='qvalue', threshold=.05 
	, conservation=0.1 , parallel=FALSE , mail=NULL , perlCommand="perl"
	,verbose=FALSE)
}
\arguments{
\item{object}{ a LowMACA class object}
\item{metric}{ a character that defines whether to use 'pvalue' or 
				'qvalue' to select significant positions. Default: 'qvalue'}
\item{threshold}{ a numeric element between 0 and 1 defining the threshold 
			of significance for the defined metric. Default: 0.05}
\item{conservation}{ a numeric value in the range of 0-1 that defines 
				the threshold of trident conservation score 
				to include the specified position. Default: 0.1}
\item{parallel}{ if TRUE, launch internal method \code{lfm} in parallel 
				mode with maximum number of available cores. If it is an integer, it sets the number of threads. 
				Default: FALSE}
\item{mail}{ if not NULL, it must be a valid email address to use EBI clustalo web service. 
			Default is to use a local clustalo installation}
\item{perlCommand}{ a character string containing the path to Perl executable. 
				if missing, "perl" will be used as default. Only used in web mode}
\item{verbose}{ logical. verbose output or not}
}
\details{
This function completes a LowMACA analysis by analyzing 
every gene or domain in the LowMACA object as a 'single sequence' analysis was started 
in the first place. The result is a dataframe showing all the significant positions 
of every gene. If you have a LowMACA object composed by 100 genes, 
it will launch 100 LowMACA single gene analyses and aggregates 
the results of every \code{lfm} launched on these 100 objects.
The output looks very similar to \code{lfm}, but in this case the 
column Multiple_Aln_pos has a different meaning. While in \code{lfm} it shows 
where the mutation falls in the consensus sequence, in this case it must be intended 
the consensus within the gene. If the original LowMACA object had mode equal to 'gene', the column 
Multiple_Aln_pos will be always equal to Amino_Acid_Position. If mode is 'pfam', it is the same unless 
a gene harbors more than one domain of the same type within its sequence. In that case, an internal alignment 
of every domain inside the protein is performed.
}
\value{
A data.frame with 10 columns corresponding to the mutations retrieved:
\enumerate{
	\item{Gene_Symbol}{ gene symbols of the analyzed genes}
	\item{Amino_Acid_Position}{ amino acidic positions relative to original protein}
	\item{Amino_Acid_Change}{ amino acid changes in hgvs format}
	\item{Sample}{ Sample barcode where the mutation was found}
	\item{Tumor_Type}{ Tumor type of the Sample}
	\item{Envelope_Start}{ start of the pfam domain in the protein}
	\item{Envelope_End}{ end of the pfam domain in the protein}
	\item{Multiple_Aln_pos}{ positions in the consensus 
						relatively to the sequence analyzed. See warnings section}
	\item{Entrez}{ entrez ids of the mutations}
	\item{Entry}{ Uniprot entry of the protein}
	\item{UNIPROT}{ other protein names for Uniprot}
	\item{Chromosome}{ cytobands of the genes}
	\item{Protein.name}{ extended protein names}
}
}
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{lfm}}
}
\examples{
#Load homeobox example
data(lmObj)
#Run lfmSingleSequence (can take some time)
if(Sys.info()[['sysname']] == 'Windows'){
	#For this example, parallelization in Windows is slower than serialization
	significant_muts <- lfmSingleSequence(lmObj , parallel=FALSE)
}else{
	significant_muts <- lfmSingleSequence(lmObj , parallel=TRUE)
}
#Show the result 
head(significant_muts)
#Show all the genes that harbor significant mutations without the alignment
unique(significant_muts$Gene_Symbol)
}
