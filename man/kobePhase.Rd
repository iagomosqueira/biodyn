\name{kobePhase}
\alias{kobePhase}
\title{kobePhase}
\arguments{
  \item{object;}{a \code{biodyn} object}
}
\value{
  A ggplot2 object
}
\description{
  produces the kobe Phase plot background, i.e. green, red
  and yellow quadrants to which layers can be added
}
\examples{
\dontrun{
    data(asp)
    kobePhase(asp)+geom_path( aes(stock,harvest)) +
    geom_point(aes(stock,harvest))
}
\seealso{
  \code{\link{kobe}}
}

