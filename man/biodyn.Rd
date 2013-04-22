\name{biodyn}
\alias{biodyn}
\title{Biomass Dynamic Model Class}
\arguments{
  \item{model}{a factor or string that specifies the model
  type, has to be one of "fox", "schaefer", "pellat",
  "gulland", "fletcher", "shepherd", "logistic", "genfit"}

  \item{params}{model parameters}
}
\value{
  biodyn object

  biodyn object
}
\description{
  A class that represents a biomass dynamic stock
  assessment model. It has slots for data and parameters
  and methods for calculating reference points etc...
  Perversely there are no methods for estimating the
  parameters from data, the class needs to be extended to
  do this , i.e. create a new class with estimation methods
  to estimate parameters

  Creates an object of the biodyn class representing a
  biomass dynamic stock assessment model.
}
\examples{
\dontrun{aspic()}
bd=biodyn("pellat",FLPar(r=0.6,k=50000,p=1,b0=1))
}

