\docType{methods}
\name{tac}
\alias{hcrParams}
\alias{tac}
\title{tac ,}
\usage{
  hcrParams(ftar, btrig, fmin, blim)
}
\arguments{
  \item{\code{object},}{an object of class \code{biodyn}
  or}

  \item{\code{harvest},}{an \code{FLQuant} object with
  harvest rate}

  \item{\code{ftar},}{an object of class \code{FLPar}}

  \item{\code{btrig},}{an object of class \code{FLPar}}

  \item{\code{fmin},}{an object of class \code{FLPar}}

  \item{\code{blim},}{an object of class \code{FLPar}}
}
\value{
  \code{FLQuant} object with TAC value(s)
}
\description{
  Calculates the Total Allowable Catch for a given harvest
  rate and stock biomass

  Combines reference points into the HCR breakpts
}
\examples{
#tac("logistic",FLPar(msy=100,k=500))
#tac("logistic",FLPar(msy=100,k=500))
}
\seealso{
  \code{\link{hcr}}, \code{\link{fwd}}

  \code{\link{hcr}}, \code{\link{hcrPlot}}
}

