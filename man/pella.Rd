\docType{methods}
\name{pella}
\alias{pella}
\title{pella}
\arguments{
  \item{\code{object},}{an object of class \code{biodyn}}

  \item{\code{index},}{an \code{FLQuant}, \code{FLQuants}
  or \code{data.frame} object with CPUE indices}

  \item{\code{cmdOps},}{a character string giving ADMB
  options see
  \url{http://www.admb-project.org/documentation/manuals/ADMBrefcard-A4.pdf/view}}
}
\description{
  Fits the Pella-Tomlinson to the catch data in a biodyn
  using CPUE index
}
\examples{
/dontrun{
data(bd)
bd=fit(bd,swonIndex)}
}

