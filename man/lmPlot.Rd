\name{lmPlot}
\alias{lmPlot}
\title{
Draw a comprehensive LowMACA plot
}
\description{
LowMACA comprehensive plot is a four layers 
plot the summarize the entire LowMACA output
}
\usage{
lmPlot(object)
}
\arguments{
\item{object}{ a LowMACA class object}
}
\details{
The method returns a plot, which is divided into four layers. 
The LowMACA object must have been passed through the methods 
\code{alignSequences} , \code{getMutations} , 
\code{mapMutations} and \code{entropy}.
The four layers of the plot are:
\enumerate{
	\item The bar plot visualized by \code{bpAll}
	\item The distribution of mutations against the 95\% confidence interval superior limit of the null hypothesis (dotted line) 
		    with orange bars representing a position with a pvalue <0.05 and a red star for qvalue<0.05
	\item The Trident score distribution
	\item The logo plot representing the consensus sequence
}
If this plot is used on a LowMACA object with a single protein, the result is formed by three layers only:
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
\code{\link{alignSequences}} \code{\link{getMutations}} \code{\link{mapMutations}} \code{\link{entropy}} \code{\link{bpAll}}
}
\examples{
#Load homeobox example and draw the plot
data(lmObj)
#Calculate statistics for nullProfile
lmObj <- entropy(lmObj)
lmPlot(lmObj)
}
