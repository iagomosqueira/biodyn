\docType{methods}
\name{bmsy}
\alias{bmsy}
\title{Biomass at Maximum Sustainable Yield.}
\arguments{
  \item{\code{object},}{an object of class \code{biodyn}
  or}

  \item{\code{object},}{a string or factor that species the
  model}

  \item{\code{params},}{an \code{FLPar} object with model
  parameters}
}
\value{
  an \code{FLPar} object with value(s) of $MSY$, $B_{MSY}$
  and $B_{MSY}$
}
\description{
  Calculates $B_{MSY}$ given the model parameters, can be
  done for a biodyn class, or by specifying the model and
  parameters
}
\examples{
refpts("logistic",FLPar(msy=100,k=500))
}
\seealso{
  \code{\link{msy}}, \code{\link{fmsy}} and
  \code{\link{bmsy}}
}

