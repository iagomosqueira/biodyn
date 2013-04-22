\docType{methods}
\name{plotSP}
\alias{plotSP}
\title{Create a \code{ggplot} plot}
\usage{
  plotSP(object,
    biomass = FLQuant(seq(0, params(object)["k"], length.out = 101)))
}
\arguments{
  \item{object,}{an object of class \code{biodyn}}

  \item{biomass,}{an object of holding biomass at beginning
  of year}
}
\value{
  an \code{ggplot2} object
}
\description{
  Creates a \code{ggplot2} object that plots time series of
  biomass, harvest rate and catch. The basic object can
  then be modified by adding ggpot2 layers.
}
\examples{
refpts("logistic",FLPar(msy=100,k=500))
}
\seealso{
  \code{\link{plotSP}}
}

