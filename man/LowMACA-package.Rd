\name{LowMACA-package}
\alias{LowMACA-package}
\alias{LowMACA}
\docType{package}
\title{
LowMACA : Low frequency Mutations Analysis via Consensus Alignment
}
\description{
The LowMACA package is a simple suite of tools to investigate and analyze 
the mutation profile of several proteins or pfam domains via consensus alignment. 
You can conduct an hypothesis driven exploratory analysis using our package 
simply providing a set of genes or pfam domains of your interest.
}
\details{
LowMACA allows to collect, align, analyze and visualize mutations from different proteins or pfam domains.
\enumerate{
\item \var{newLowMACA:} construct a LowMACA object with your proteins or pfam
\item \var{setup:} align sequences, get mutations and map mutations on the consensus sequence
\item \var{entropy:} calculate entropy score and pvalues for every position
\item \var{lfm:} retrieve significant position
\item \var{lmPlot:} visualize mutations on the consensus sequence, conservation and significant clusters
}
}
\author{
Stefano de Pretis , Giorgio Melloni

Maintainer: <ste.depo@gmail.com> <melloni.giorgio@gmail.com>
}
\references{
doi:10.1186/gm563 923 Melloni et al.: 
DOTS-Finder: a comprehensive tool for assessing driver genes in cancer genomes. 
Genome Medicine 2014 6:44
}
\keyword{ package }
\seealso{
\href{https://cgsb.genomics.iit.it/wiki/projects/LowMACA}{LowMACA project website}
}
\examples{
#Create an object of class LowMACA for RAS domain family
lm <- newLowMACA(pfam="PF00071" , genes=c("KRAS" , "NRAS" , "HRAS"))
#Select melanoma, breast cancer and colorectal cancer
lmParams(lm)$tumor_type <- c("skcm" , "brca" , "coadread")
#Align sequences, get mutation data and map them on consensus
lm <- setup(lm)
#Calculate statistics
lm <- entropy(lm)
#Retrieve original mutations
lfm(lm)
#Plot
bpAll(lm)
lmPlot(lm)
protter(lm)
}
