\name{newLowMACA}
\alias{newLowMACA}
\title{
Construct a LowMACA object
}
\description{
Constructor for the class LowMACA. It initializes a LowMACA object with default parameters
}
\usage{
newLowMACA(genes = NULL, pfam = NULL)
}
\arguments{
\item{genes}{ a character vector of gene symbols in HGNC format or a integer vector of Entrez IDs.
          If \code{pfam} is defined, it can be set to NULL}
\item{pfam}{ a character vector of pfam IDs.
           If \code{genes} is defined, it can be set to NULL}
}
\details{
When a LowMACA object is initialized, the \code{arguments} slot is filled 
with the input data and default parameters and path to clustalomega aligner. 
See \code{lmParams} and \code{parallelize} to change them.
}
\value{
An object of class "LowMACA". See \code{\link{LowMACA-class}}
}
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{lmParams}} \code{\link{parallelize}}
}
\examples{
#Set Genes and pfam for the analysis
Genes <- c("ADNP","ALX1","ALX4","ARGFX","CDX4","CRX"
    	,"CUX1","CUX2","DBX2","DLX5","DMBX1","DRGX"
			,"DUXA","ESX1","EVX2","HDX","HLX","HNF1A"
			,"HOXA1","HOXA2","HOXA3","HOXA5","HOXB1","HOXB3"
			,"HOXD3","ISL1","ISX","LHX8")
Pfam <- "PF00046"
#LowMACA object of pfam PF00046 filtered by Genes
lm <- newLowMACA(genes=Genes , pfam=Pfam)
#LowMACA object of the entire pfam PF00046
lm <- newLowMACA(pfam=Pfam)
#LowMACA object of entire canonical proteins associated to Genes
lm <- newLowMACA(genes=Genes)
}
