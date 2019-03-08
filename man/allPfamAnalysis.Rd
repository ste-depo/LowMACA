\name{allPfamAnalysis}
\alias{allPfamAnalysis}
\title{
Global analysis of a repository of mutations
}
\description{
Given a repository of mutations, the method \code{allPfamAnalysis} 
launches the analysis of all the Pfams and single sequences which 
are involved with at least one mutation.
}
\usage{
allPfamAnalysis(repos
, allLowMACAObjects=NULL
, mutation_type=c("missense", "all", "truncating" , "silent") 
, NoSilent=TRUE
, mail=NULL
, perlCommand="perl"
, verbose=FALSE
, conservation=0.1
, use_hmm=FALSE
, datum=FALSE
, clustal_cmd="clustalo"
, BPPARAM=bpparam("SerialParam"))
}
\arguments{
\item{repos}{ either a data.frame or a filename containing the data to analyze}
\item{allLowMACAObjects}{ filename of a RData file to save all 
							the LowMACA object \code{allPfamsLM} produced by the function. It can be usefull for plotting
							a specific Pfam after the analysis, but it can be a pretty large object. Default NULL}
\item{mutation_type}{ type of mutation to be considered for the analysis. 
		Default to missense.}
\item{NoSilent}{ logical indicating if Silent mutations should be deleted or not.
	Default TRUE}
\item{mail}{ if not NULL, it must be a valid email address to use EBI 
	clustalo web service. Default is to use a local clustalo installation}
\item{perlCommand}{ a character string containing the path to Perl executable. 
				if missing, "perl" will be used as default. Only used if mail is set}
\item{verbose}{ logical. verbose output or not}
\item{conservation}{ a number between 0 and 1. Represents the minimum level of conservation to test a mutation}
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
\item{clustal_cmd}{ path to clustalomega executable. default is to check "clustalo" in the PATH}
\item{BPPARAM}{ An object of class \code{\link[BiocParallel]{BiocParallelParam-class}} specifiying parameters related to
  the parallel execution of some of the tasks and calculations within this function.
  See function \code{\link[BiocParallel]{register}} from the \code{BiocParallel} package.}
}
\details{
This function takes a data.frame or a tab delimited text file in LowMACA format (see \code{\link{LowMACA_AML}})
and perform a full analysis of the dataset. It basically divide the mutations into their Pfam and launch many LowMACA 
analysis as many Pfam are hit by mutations up to the \code{lfm} function. Every significant position after \code{lfm} 
is tested at gene level. A binomial test is performed to see if the ratio between the number of mutations 
in the significant position over the total number of mutations is higher than expected by chance at gene level.
The significant mutations of all the \code{lfm} functions are aggregated in one single data.frame.
}
\value{
A list of two dataframes named 'AlignedSequence' and 'SingleSequence'

The first dataframe is the result of the alignment based analysis. 
Every gene is aggregated by its corresponding Pfam domain.
	\item{Gene_Symbol}{ gene symbols of the analyzed genes}
	\item{Multiple_Aln_pos}{ positions in the consensus 
						relatively to the sequence analyzed.}
	\item{Pfam_ID}{ Pfam name analyzed}
	\item{binomialPvalue}{ pvalue of the single gene test, See details}
	\item{Amino_Acid_Position}{ amino acidic positions relative to original protein}
	\item{Amino_Acid_Change}{ amino acid changes in hgvs format}
	\item{Sample}{ Sample barcode where the mutation was found}
	\item{Tumor_Type}{ Tumor type of the Sample}
	\item{Envelope_Start}{ start of the pfam domain in the protein}
	\item{Envelope_End}{ end of the pfam domain in the protein}
	\item{metric}{ qvalue of the position in the multiple alignment of Pfam domains}
	\item{Entrez}{ entrez ids of the mutations}
	\item{Entry}{ Uniprot entry of the protein}
	\item{UNIPROT}{ other protein names for Uniprot}
	\item{Chromosome}{ cytobands of the genes}
	\item{Protein.name}{ extended protein names}
The second dataframe represent the result of LowMACA on every couple gene-domain 
when it is not aligned with any other member of the same Pfam ID.
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
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{lfm}}, \code{\link{LowMACA_AML}}
}
\examples{
#Load Homeobox example
data(lmObj)
#Extract the data inside the object as a toy example
myData <- lmMutations(lmObj)$data
#Run allPfamAnalysis on every mutations
significant_muts <- allPfamAnalysis(repos=myData)
#Show the result of alignment based analysis
head(significant_muts$AlignedSequence)
#Show all the genes that harbor significant mutations
unique(significant_muts$AlignedSequence$Gene_Symbol)
#Show the result of the Single Gene based analysis
head(significant_muts$SingleSequence)
#Show all the genes that harbor significant mutations
unique(significant_muts$SingleSequence$Gene_Symbol)
}
