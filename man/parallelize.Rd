\name{parallelize}
\alias{parallelize}
\alias{parallelize<-}
\title{
Show and set parallelization options
}
\description{
Method for objects of class LowMACA. It can show parallelization parameters 
of an object of class LowMACA and switch off and on parallelization 
of \code{\link{alignSequences}} and \code{\link{getMutations}} method
}
\usage{
parallelize(object)
parallelize(object) <- value
}
\arguments{
\item{object}{ object of class LowMaca}
\item{value}{ a named list containing logical values.
              Default list(getMutations=FALSE , makeAlignment=TRUE)}
}
\details{
With \code{getMutations=TRUE}, the \code{\link{getMutations}} 
method runs in parallel during the queries to the different tumor_types. 
This can result in an overload to the cBioPortal database and the function returns error.
With \code{makeAlignment=TRUE}, clustalo should run in parallel. 
Nevertheless, clustalo can be parallelized only if the OpenMP C library is correctly functioning. 
}
\value{
If \code{parallelize} is used as a show method 
it returns a named list of two elements: \code{getMutations}
and \code{makeAlignment}
}
\author{
Stefano de Pretis , Giorgio Melloni
}
\seealso{
\code{\link{getMutations}}
}
\examples{
#Construct a LowMACA object
lm <- newLowMACA(pfam="PF12906")
#Show parallelize default
parallelize(lm)
#Change all parameters
parallelize(lm) <- list(getMutations=TRUE , makeAlignment=FALSE)
#Change just one parameter
parallelize(lm)[['getMUtations']] <- TRUE
}
