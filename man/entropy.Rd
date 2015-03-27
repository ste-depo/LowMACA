\name{entropy}
\alias{entropy}
\title{
Calculate LowMACA statistics
}
\description{
\code{entropy} is a method for objects of class LowMACA. 
It calculates global entropy score of the mutation profile of the alignment 
and a test for every position in the consensus comparing the number 
of observed mutations against a weigthed random uniform distribution.
}
\usage{
entropy(object, bw = NULL)
}
\arguments{
\item{object}{ an object of class LowMACA}
\item{bw}{ a character string or a numeric positive value representing the desired bandwith 
			to launch the function density for the uniform distribution. 
			0 will not launch density (every position is not aggregated to the surrounded ones) , 
			'auto' will let the simulation decide according to the Silverman's rule of thumb 
			and every other number is a user defined bandwidth passed to the function \code{density}.}
}
\value{
\code{entropy} returns an object of class LowMACA 
updating the slot \code{entropy} and the slot \code{alignment} 
in the \code{df} data.frame by adding 6 new columns
\describe{
alignment slot update:
\itemize{
	\item{mean}{ a numeric vector representing the mean value of the empirical 
		uniform function at every position in the consensus}
	\item{lTsh}{ a numeric vector representing the limit inferior of the 95\% confidence interval of the 
		empirical uniform function at every position in the consensus}
	\item{uTsh}{ a numeric vector representing the limit superior of the 95\% confidence interval of the 
		empirical uniform function at every position in the consensus}
	\item{profile}{ a numeric vector representing the density of mutations at every position in the sample 
			normalized by the number of position. In case of bandwidth 0, this vector is equal 
			to the number of mutations divided by the total number of mutations}
	\item{pvalue}{ a numeric vector representing the pvalue of the number of mutations 
			found at every position against the weigthed random uniform distribution of mutations}
	\item{qvalue}{ a numeric vector representing the corrected pvalues using FDR method. 
		Only positions with a conservation score >= 10\% are considered}
}
}
}
\details{
The parameter \code{bw} overwrites the bandwidth set with \code{lmParams}. Therefore, 
if bw is set to NULL, the method entropy uses the predefined bandwidth of the LowMACA object.
}
\references{
doi:10.1186/gm563 923 Melloni et al.: 
\emph{DOTS-Finder: a comprehensive tool for assessing driver genes in cancer genomes.} 
Genome Medicine 2014 6:44

Silverman, B. W. (1986) \emph{Density Estimation}. London: Chapman and Hall.
}
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{alignSequences}} \code{\link{lmParams}} \code{\link{lmEntropy}}
}
\examples{
#Load homeobox example and run entropy
data(lmObj)
lmObj <- entropy(lmObj)
lmEntropy(lmObj)
}
