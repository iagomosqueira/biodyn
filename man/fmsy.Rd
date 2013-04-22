\docType{methods}
\name{fmsy}
\alias{fmsy}
\title{Fishing Mortality at Maximum Sustainable Yield.}
\arguments{
  \item{\code{object},}{an object of class \code{biodyn}
  or}

  \item{\code{object},}{a string or factor that species the
  model}

  \item{\code{params},}{an \code{FLPar} object with model
  parameters}
}
\value{
  an \code{FLPar} object with value(s) of $F_{MSY}$
}
\description{
  Calculates $F_{MSY}$ given the model parameters, can be
  done for a biodyn class, or by specifying the model and
  parameters
}
\examples{
fmsy("logistic",FLPar(msy=100,k=500))
}
\seealso{
  \code{\link{msy}}, \code{\link{bmsy}} and
  \code{\link{refpts}}
}

