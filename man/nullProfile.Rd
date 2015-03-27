\name{nullProfile}
\alias{nullProfile}
\title{
Draw a mutational profile plot
}
\description{
\code{nullProfile} is a method for objects of class LowMACA 
that draw a barplot highlighting the significant clusters 
of mutations found by LowMACA
}
\usage{
nullProfile(object)
}
\arguments{
\item{object}{ an object of class LowMACA}
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
