\docType{methods}
\name{hcr}
\alias{hcr}
\title{HCR}
\arguments{
  \item{\code{object}}{an object of class \code{biodyn} or}

  \item{\code{params}}{\code{FLPar} object with hockey
  stick HCR parameters}

  \item{\code{yrs},}{numeric vector with yrs for HCR
  prediction}

  \item{\code{refYrs},}{numeric vector with years used to
  for stock/ssb in HCR}

  \item{\code{tac},}{\code{logical} should return value be
  TAC rather than F?}

  \item{\code{bndF},}{\code{vector} with bounds on
  iter-annual variability on F}

  \item{\code{bndTac},}{\code{vector} with bounds on
  iter-annual variability on TAC}
}
\value{
  \code{FLPar} object with value(s) for F or TAC if
  tac==TRUE
}
\description{
  Harvest Control Rule, calculates F, or Total Allowable
  Catch (TAC) based on a hockey stock harvest control rule.
}
\examples{
#hcr("logistic",FLPar(msy=100,k=500))
}
\seealso{
  \code{\link{bmsy}}, \code{\link{fmsy}},
  \code{\link{fwd}}, \code{\link{hcr}} and \code{\link{hcrParams}}
}

