\name{lmEntropy}
\alias{lmEntropy}
\title{
Show Entropy Information Contained in a LowMACA object
}
\description{
Method for objects of class LowMACA. It can show the results of entropy analysis
performed on the LowMACA object by the function \code{\link{entropy}}
}
\usage{
lmEntropy(object)
}
\arguments{
\item{object}{object of class LowMaca}
}
\value{
A list containing the following elements:
\itemize{
  \item{bw}{ a numeric value that represents the bandwidth used to calculate the Shannon entropy score}
  \item{uniform}{ an object of class function that was used to calculate the score}
  \item{absval}{ a numeric value representing the Shannon entropy of the sample data}
  \item{log10pval}{ a numeric value representing the pvalue of the Shannon entropy score against 
        a gamma distribution with same mean and variance as the empirical 
        uniform distribution in -log10 scale}
  \item{pvalue}{ a numeric value representing the pvalue of the Shannon entropy score against 
      a gamma distribution with same mean and variance as the empirical 
      uniform distribution}
}
}
\author{
Stefano de Pretis, Giorgio Melloni
}
\seealso{
\code{\link{entropy}}
}
\examples{
data('lmObj')
lmObj <- entropy(lmObj)
lmEntropy(lmObj)
}
