\name{calcK}
\alias{calcK}
\title{Calculate Carrying Capacity}
\usage{
  calcK(msy, r = 0.6, p = 1, params = FLPar(r = r, p = p))
}
\arguments{
  \item{\code{msy},}{a guess for $MSY$}

  \item{\code{r},}{a guess for $r$ the population growth
  rate}

  \item{\code{p},}{a guess for $p$ the shape parameter}
}
\value{
  an \code{FLPar} object with an estimate for $k$
}
\description{
  Calculates $k$ given msy, r and K for a Pella-Tomlinson
  biomass dynamic model
}
\examples{
calcK(5000,r=.6,p=1.0)
}
\seealso{
  \code{\link{msy}} and \code{\link{bmsy}}
}

