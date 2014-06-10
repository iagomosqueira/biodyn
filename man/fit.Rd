\docType{methods}
\name{fit}
\alias{fit}
\alias{profile}
\title{fit}
\arguments{
  \item{\code{object},}{an object of class \code{biodyn}}

  \item{\code{index},}{an \code{FLQuant}, \code{FLQuants}
  or \code{data.frame} object with CPUE indices}

  \item{\code{cmdOps},}{a character string giving ADMB
  options see
  \url{http://www.admb-project.org/documentation/manuals/ADMBrefcard-A4.pdf/view}}

  \item{\code{fitted},}{an object of class \code{biodyn}}

  \item{\code{index},}{an \code{FLQuant}, \code{FLQuants}
  or \code{data.frame} object with CPUE indices}

  \item{\code{cmdOps},}{a character string giving ADMB
  options see
  \url{http://www.admb-project.org/documentation/manuals/ADMBrefcard-A4.pdf/view}}
}
\description{
  Estimates parameters in a \code{biodyn} class by fitting
  catch to CPUE indeices

  Profiles biodyn
}
\examples{
/dontrun{
data(bd)
bd=fit(bd,swonIndex)}
/dontrun{
data(asp)
cpue=FLQuant(asp
}

