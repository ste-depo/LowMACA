\name{lfm}
\alias{lfm}
\title{
Show significant clusters of mutations
}
\description{
The method \code{lfm} (low frequency mutations) retrieve 
the original mutations that created the significant clusters 
calculated with \code{entropy} on the consensus
}
\usage{
lfm(object , metric='qvalue', threshold=.05, conservation=NULL)
}
\arguments{
\item{object}{ a LowMACA class object}
\item{metric}{ a character that defines whether to use 'pvalue' or 
				'qvalue' to select significant positions. Default: 'qvalue'}
\item{threshold}{ a numeric defining the threshold of significance for the defined metric. Default: 0.05}
\item{conservation}{ a numeric value in the range of 0-1 that defines 
				the threshold of trident conservation score to include the specified position.
				The default value is inherited from the slot entropy, whose default is 0.1}
}
\details{
After the alignment, we lose every information about the original sequences used as input.
The consensus sequence is in fact an alignment that could not represent the reality of human proteins.
\code{lfm} allows to go back on the original dataset and retrieve the proteins and the real positions
of the mutations that we consider 'conserved'.
}
\value{
A data.frame with 13 columns corresponding to the mutations retrieved:
\enumerate{
	\item{Gene_Symbol}{ gene symbols of the mutations}
	\item{Amino_Acid_Position}{ amino acidic positions relative to original protein}
	\item{Amino_Acid_Change}{ amino acid changes in hgvs format}
	\item{Sample}{ Sample barcode where the mutation was found}
	\item{Tumor_Type}{ Tumor type of the Sample}
	\item{Envelope_Start}{ start of the pfam domain in the protein}
	\item{Envelope_End}{ end of the pfam domain in the protein}
	\item{Multiple_Aln_pos}{ positions in the consensus}
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
\code{\link{entropy}}
}
\examples{
#Load homeobox example and launch entropy method
data(lmObj)
lmObj <- entropy(lmObj)
significant_muts <- lfm(lmObj)
#Display original mutations that formed significant clusters (column Multiple_Aln_pos)
head(significant_muts)
#Position 4 has a qvalue<0.05
#What are the genes mutated in position 4 in the consensus?
cluster_4_genes <- significant_muts[ significant_muts[['Multiple_Aln_pos']]==4 , 'Gene_Symbol']
#Display the genes and their number of mutation in consensus position 4
sort(table(cluster_4_genes))
}
