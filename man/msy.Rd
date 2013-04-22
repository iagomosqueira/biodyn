\docType{methods}
\name{msy}
\alias{msy}
\title{Maximum Sustainable Yield.}
\arguments{
  \item{\code{object},}{an object of class \code{biodyn}
  or}

  \item{\code{object},}{a string or factor that species the
  model}

  \item{\code{params},}{an \code{FLPar} object with model
  parameters}
}
\value{
  an \code{FLPar} object with value(s) of MSY
}
\description{
  Calculates MSY given the model parameters, can be done
  for a biodyn class, or by specifying the model and
  parameters
}
\examples{
\dontrun{ msy("logistic",FLPar(msy=100,k=500))}
}
\seealso{
  \code{\link{bmsy}}, \code{\link{fmsy}} and
  \code{\link{refpts}}
}

