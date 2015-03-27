\name{bpAll}
\alias{bpAll}
\title{
Draw a mutation barplot
}
\description{
bpAll draws a stacked barplot of the mutations mapped on the consensus sequence
}
\usage{
bpAll(object)
}
\arguments{
  \item{object}{ an object of class LowMACA}
}
\details{
Returns a barplot in which mutations are stacked per position
on the consensus sequence. Every color represent mutations taht map
on the same input sequence (either a protein or a pfam)
The LowMACA object must pass through the methods 
\code{alignSequences} , \code{getMutations} , \code{mapMutations}
}
\value{
NULL
}
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{lmPlot}}
}
\examples{
#Load homeobox example and draw plot
data(lmObj)
lmObj <- entropy(lmObj)
bpAll(lmObj)
}
