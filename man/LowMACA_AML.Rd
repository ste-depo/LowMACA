\name{LowMACA_AML}
\alias{LowMACA_AML}
\docType{data}
\title{
Example of a LowMACA object
}
\description{
A data frame containing TCGA AML data in the format accepted by LowMACA
}
\usage{data("LowMACA_AML")}
\format{
A data.frame of 8 columns:
\enumerate{
	\item{Entrez}{ gene ID number}
	\item{Gene_Symbol}{ HGNC official gene symbol}
	\item{Amino_Acid_Letter}{ original amino acid letter in the 
							position of the mutation}
	\item{Amino_Acid_Position}{ position of the mutation 
								relative to the protein}
	\item{Amino_Acid_Change}{ amino acid change in hgvs format, like G12V}
	\item{Mutation_Type}{ classification of mutation according to MAF format.}
	\item{Sample}{ name of the sample where the mutation was found}
	\item{Tumor_Type}{ type of tumor, if applicable}
}
}
\source{
Adapted from \href{https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/laml/}{TCGA ftp repository}
}
\seealso{
\href{https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format}{MAF format specification}
\href{http://www.hgvs.org/mutnomen/}{HGVS}
}
\examples{
#Load LowMACA_AML and show its structure
data(LowMACA_AML)
str(LowMACA_AML)
}
\keyword{datasets}
