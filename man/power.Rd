\docType{methods}
\name{power}
\alias{power}
\title{power}
\arguments{
  \item{\code{object},}{an object of class \code{biodyn}}

  \item{\code{test},}{a \code{data.frame} with 3 columns
  \code{param,q,lower.tail} for the choice of parameter,
  the desired qunatile and teh tail of the distribution for
  the test.}
}
\value{
  an \code{FLPar} object with probability of being greater
  (lower.tail=FALSE) or less (lower.tail=TRUE) then q
}
\description{
  Estimates parameters in a \code{biodyn} class by fitting
  catch to CPUE indeices
}

