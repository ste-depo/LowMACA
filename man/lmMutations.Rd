\name{lmMutations}
\alias{lmMutations}
\title{
Show Mutation Data Contained in a LowMACA object
}
\description{
Method for objects of class LowMACA. It can show the mutation data contained 
within the LowMACA object that has been retrieved from \code{\link{getMutations}}
method.
}
\usage{
lmMutations(object)
}
\arguments{
\item{object}{object of class LowMaca}
}
\value{
A list containing the following elements:
\itemize{
  \item{data}{ a data.frame describing the mutations 
            on every genes and their effect the amino acids they belong to}
  \item{freq}{ a data.frame containing the absolute number of mutated patients by gene and selected tumor types 
              (this is useful to explore the mutational landscape of your genes in the different tumor types)}
  \item{aligned}{ a matrix where rows represent proteins/pfam, and columns report 
                the number of mutations on every position of the  consensus}
}
}
\author{
Stefano de Pretis, Giorgio Melloni
}
\seealso{
\code{\link{getMutations}}
}
\examples{
data('lmObj')
str(lmMutations(lmObj))
}
