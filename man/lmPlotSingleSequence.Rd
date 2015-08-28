\name{lmPlotSingleSequence}
\alias{lmPlotSingleSequence}
\title{
Draw a LowMACA comprehensive plot of a 
specified gene within a LowMACA object
}
\description{
LowMACA comprehensive plot is a four layers 
plot that summarize the entire LowMACA output
}
\usage{
lmPlotSingleSequence(object , gene , mail=NULL , perlCommand="perl")
}
\arguments{
\item{object}{ a LowMACA class object}
\item{gene}{ a Gene Symbol that identifies 
		one of the gene analyzed in the LowMACA object}
\item{mail}{ if not NULL, it must be a valid email address to use EBI clustalo web service. 
			Default is to use a local clustalo installation}
\item{perlCommand}{ a character string containing the path to Perl executable. 
				if missing, "perl" will be used as default. Only used in web mode}
}
\details{
If the specified gene has more than one domain of the same type and mode is pfam, 
the plot is composed by four layers:
\enumerate{
	\item The bar plot visualized by \code{bpAll}
	\item The distribution of mutations against the 95\% confidence interval superior limit of the null hypothesis (dotted line) 
		    with orange bars representing a position with a pvalue <0.05 and a red star for qvalue<0.05
	\item The Trident score distribution
	\item The logo plot representing the consensus sequence
}
If the specified gene has only one domain of the same type and mode is pfam, 
the plot is composed by two layers:
\enumerate{
	\item The bar plot visualized by \code{bpAll}
	\item The distribution of mutations against the 95\% confidence interval superior limit of the null hypothesis (dotted line) 
		    with orange bars representing a position with a pvalue <0.05 and a red star for qvalue<0.05
}
If mode is gene, the plot is composed by three layers:
\enumerate{
	\item The bar plot visualized by \code{bpAll}
	\item The Pfam domains structure inside the protein
	\item The distribution of mutations against the 95\% confidence interval superior limit of the null hypothesis (dotted line) 
		    with orange bars representing a position with a pvalue <0.05 and a red star for qvalue<0.05
}
}
\value{
NULL
}
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{lmPlot}} \code{\link{bpAll}}
}
\examples{
#Load homeobox example and draw the plot
data(lmObj)
#DUXA has a significant cluster of mutation
#Plot Mutations on DUXA gene in the 
#original sequences of its domains PF00046
lmPlotSingleSequence(lmObj , gene="DUXA")
}
