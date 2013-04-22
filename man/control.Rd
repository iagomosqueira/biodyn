\docType{methods}
\name{setControl<-}
\alias{setControl<-}
\title{control}
\arguments{
  \item{\code{object},}{an object of class \code{biodyn}}
}
\description{
  Sets up the control slot in a biodyn object given the
  values in the \code{params} slot. The starting values are
  set to the values in \code{params} and the min and max
  bounds to .1 and 10 times of these.
}
\examples{
/dintrun{
data(bd)
control(bd) <-params(bd)
params(bd)
control(bd)}
}

