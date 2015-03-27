\name{getMutations}
\alias{getMutations}
\title{
Retrieve mutation data for a LowMACA object
}
\description{
Exploting the capabilities of the cgdsr package, this method downloads and parse 
the mutation data of the specified genes in the selected tumor types. 
It also aggregates and show the frequencies of mutations of every gene in the different tumor types.
}
\usage{
getMutations(object, repos = NULL)
}
\arguments{
\item{object}{ a LowMACA class object}
\item{repos}{ a data.frame containing mutations for the specified genes 
                in the LowMACA object in case of custom mutation data. Default NULL}
}
\details{
With \code{repos=NULL}, the method is a wrapper around \code{\link[cgdsr]{getMutationData}} 
method from package \code{\link[cgdsr]{cgdsr-package}}. The output of the method is moduled 
by the parameters in \code{lmParams("LowMACA_object")}. See \code{\link{lmParams}} for further information.
}
\value{
An object of class LowMACA is returned with an update in the slot \code{mutations}.
See \code{\link{lmMutations}} method.
}
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{lmParams}} \code{\link[cgdsr]{getMutationData}} \code{\link{lmMutations}}
}
\examples{
#Create an object of class LowMACA
lm <- newLowMACA(pfam="PF12906")
#Change some paramters
#By default, LowMACA retrieve only missense mutations.
#We want all mutations
lmParams(lm)[['mutation_type']] <- 'all'
#By default, LowMACA takes mutations from all the kinds of tumor
#We want just prostate cancer samples
lmParams(lm)[['tumor_type']] <- 'prad'
lm <- getMutations(lm)
}
