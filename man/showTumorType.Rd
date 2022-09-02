\name{showTumorType}
\alias{showTumorType}
\title{
List of tumor types
}
\description{
Show all the possible tumor types accepted by LowMACA
}
\usage{
showTumorType()
}
\value{
A named vector of all the tumor types available in cBioportalData package 
that can be passed to the method lmParams. Every element is the 
aggregation of all the available sequenced data from all the studies 
involved in a particular tumor type.
}
\details{
This method is a wrapper around \code{\link[cBioportalData]{cBioportalData-getStudies}} 
and show all the barcodes for the tumor types as used by cBioPortal.
}
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{lmParams}} \code{\link[cBioportalData]{cBioportalData-getStudies}}
}
\examples{
data('lmObj')
out <- showTumorType()
chosenTumors <- out[1:3]
lmParams(lmObj)$tumor_type <- chosenTumors
}
