\name{mapMutations}
\alias{mapMutations}
\title{
Map mutations on consensus sequence
}
\description{
\code{mapMutations} is a method for the class LowMACA 
that re-maps the mutations on a sequence to the relative position in a consensus sequence.
}
\usage{
mapMutations(object)
}
\arguments{
\item{object}{ an object of class LowMACA}
}
\details{
Every position in the consensus alignement correspond to different positions 
in the single aligned sequences. The mutations are mapped according to this 
scheme that can be evinced from the slot \code{alignment}. 
\code{mapMutations} must be called after \code{alignSequences} and \code{getMutations}
}
\value{
An object of class LowMACA with an update in the slot \code{mutations}. 
\code{mapMutations} add a object named \code{aligned} of class matrix in this slot 
that represents the absolute number of mutations in each sequence/position 
in the consensus as a matrix.
}
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{getMutations}} \code{\link{alignSequences}} \code{\link{LowMACA-class}}
}
\examples{
#Create an object of class LowMACA
lm <- newLowMACA(pfam="PF12906")
#Align the sequences, requires clustalo
\dontrun{lm <- alignSequences(lm)}
#Get mutations from the corresponding genes
\dontrun{lm <- getMutations(lm)}
#Map mutations on the consensus sequence
\dontrun{lm <- mapMutations(lm)}
}
