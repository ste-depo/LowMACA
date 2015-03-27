\name{protter}
\alias{protter}
\title{
Draw a Protter plot
}
\description{
This is a wrapper around Protter web service for LowMACA class objects that draw a protter style plot.
}
\usage{
protter(object, filename = "protter.png", threshold = 0.05)
}
\arguments{
\item{object}{ an object of class LowMACA}
\item{filename}{ a character string that identifies the file name
 				where protter plot will be stored. Default "protter.png"}
\item{threshold}{ a numeric value in the interval (0 , 1] 
				that identifies the significant mutations. Default 0.05}
}
\details{
Using the information in the slot \code{alignment}, a request is send to Protter server. 
Protter will predict a possible sencondary structure for the consensus sequence 
(if possible) and highlights the significant clusters of mutations found by LowMACA (if any). 
A significant pvalue is colored in orange, a significant qvalue is colored in red.
}
\value{
NULL
}
\references{
\href{http://wlab.ethz.ch/protter/start/}{Protter website}
}
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{LowMACA-class}} \code{\link{entropy}}
}
\examples{
#Load homeobox example
data(lmObj)
#Calculate statistics
lmObj <- entropy(lmObj)
#Create protter.png
protter(lmObj)
}
