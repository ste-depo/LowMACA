\name{lmAlignment}
\alias{lmAlignment}
\title{
Show Alignment Results from a LowMACA object
}
\description{
Method for objects of class LowMACA. It can show the results of the alignment
procedure that has been performed on the LowMACA object
}
\usage{
lmAlignment(object)
}
\arguments{
\item{object}{object of class LowMaca}
}
\value{
A list containing the following elements:
\itemize{
  \item{ALIGNMENT}{ an object of class data.frame containing the mapping 
                  of the position of the original amino acids to the consensus sequence}
  \item{SCORE}{ a list of two objects
    \itemize{
      \item{DIST_MAT}{ a matrix of the pairwise similarities between sequences 
                      as resulted after the multiple alignment (from 0\% to 100\%)}
      \item{SUMMARY_SCORE}{ a data.frame containing summary descriptives of the distance matrix}
      \item{CLUSTAL}{ an object of class "AAMultipleAlignment" as provided by Biostrings R package}
      \item{df}{ a dataframe containing the predicted consesus sequence and the trident conservation score at every position}
    }
  }
}
}
\author{
Stefano de Pretis, Giorgio Melloni
}
\seealso{
\code{\link{alignSequences}}
}
\examples{
data('lmObj')
str(lmAlignment(lmObj))
}
