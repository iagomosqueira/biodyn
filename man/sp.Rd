\docType{methods}
\name{computeSP}
\alias{computeSP}
\title{Calculates surplus production}
\arguments{
  \item{\code{object},}{an object of class \code{biodyn}}

  \item{\code{biomass},}{stock biomaas, may be a
  \code{numerix}, \code{FLQuant} or missing. In the latte
  case the stock slot will be used.}
}
\value{
  an \code{FLPar} object
}
\description{
  Calculates the surplus production for a biomass dynamic
  model given a level of stock biomass
}
\examples{
\dontrun{ computeSP(bd,seq(0,params(bd)["k"])) }
}
\seealso{
  \code{\link{plotSP}}
}

