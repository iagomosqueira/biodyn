\docType{methods}
\name{refpts}
\alias{refpts}
\title{Maximum Sustainable Yield reference points}
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
  Calculates $MSY$, $B_{MSY}$ and $F_{MSY}$ given the model
  parameters, can be done for a biodyn class, or by
  specifying the model and parameters
}
\examples{
refpts("logistic",FLPar(msy=100,k=500))
}
\seealso{
  \code{\link{msy}}, \code{\link{bmsy}} and
  \code{\link{refpts}}
}

