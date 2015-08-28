\name{BLOSUM62}
\alias{BLOSUM62}
\docType{data}
\title{
BLOSUM62 matrix
}
\description{
A substitution matrix used for sequence alignment of proteins.
In LowMACA, it is used to calculate the trident conservation score.
}
\usage{data("BLOSUM62")}
\format{
A squared numeric matrix with aminoacids as rownames and colnames
}
\source{
\href{http://www.uky.edu/Classes/BIO/520/BIO520WWW/blosum62.htm}{BLOSUM62 from NCBI}
}
\examples{
#Load BLOSUM62 and show its structure
data(BLOSUM62)
str(BLOSUM62)
}
\keyword{datasets}
