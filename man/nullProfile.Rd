\name{nullProfile}
\alias{nullProfile}
\title{
Draw a mutational profile plot
}
\description{
\code{nullProfile} is a method for objects of class LowMACA 
that draw a barplot highlighting the significant clusters 
of mutations found by LowMACA statistics
}
\usage{
nullProfile(object , conservation=NULL, windowlimits=NULL)
}
\arguments{
\item{object}{ an object of class LowMACA}
\item{conservation}{ a numeric value in the range of 0-1 that defines 
				the threshold of trident conservation score to include the specified position.
				The default value is inherited from the slot entropy, whose default is 0.1}
\item{windowlimits}{ A vector indicating which amino acids residues will be 
				plotted. The vector refers to the positions in the global
				alignment. By default this parameter is set to NULL, 
				that means that all the amino acids will be displayed. }
}
\details{
This method draw the second layer of the \code{lmPlot} of LowMACA. 
The blue dotted line is a curve that pass through all the points of 
the upper limit of the 95\% confidence interval for the single position 
test performed by \code{entropy} (one point per position in the consensus). 
The black bars represents the density of mutations in our sample. 
If a bar passes the blue line, it will be depicted in orange (significant pvalue). 
After the correction for multiple testing, red stars appears at the top of the orange bars 
if a cluster is below 0.05 for the qvalue and has a conservation trident score of at least 0.1.
}
\value{
NULL
}
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{lmPlot}} \code{\link{entropy}}
}
\examples{
#Load homeobox example
data(lmObj)
#Calculate statistics
lmObj <- entropy(lmObj)
nullProfile(lmObj)
}
