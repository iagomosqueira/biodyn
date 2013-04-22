\docType{methods}
\name{fwd}
\alias{fwd}
\title{Forward Projections}
\arguments{
  \item{\code{object},}{an object of class \code{biodyn}
  or}

  \item{\code{ctrl},}{an \code{FLQuants} object with names
  corresponding to "f", "catch" etc.}
}
\description{
  Projects the biomass dynamic model into the future for a
  given catch, harvest rate etc.
}
\examples{
bd=biodyn("logistic",FLPar(msy=1000,k=50000,b0=1),stock=FLQuant(50000,dimnames=list(year=1:50)))
 bd.=fwd(bd,catch=FLQuant(500,dimnames=list(year=1:50)))
}

