\name{biodyn-package}
\alias{biodyn-package}
\docType{package}
\title{
Biomass Dynamic Model
}
\description{
Implements a biomass dynamic model as an R class for fish stock abundance. The class can be used to fit catch and catch per unit effort to estimate stock abundance. It also has methods for evaluating the goodness of fit, estimating reference points and for modelling harvest control rules.

}
\details{
\tabular{ll}{
Package: \tab biodyn\cr
Version: \tab 0.0.1\cr
Depends: \tab 
R(>= 2.15.0),
ggplot2,
kobe,
diags,
FLCore(>= 2.2),
FLash,
FLAdvice,
ggplotFL\cr
Collate: \tab 
'biodyn-class.R'
'biodyn-generic.R'
'biodyn-control.R'
'biodyn-accessors.R'
'biodyn-constructors.R'
'biodyn-OEM.R'
'biodyn-admb.R'
'biodyn-coerce.R'
'biodyn-methods.R'
'biodyn-diags.R'
'biodyn-fwd.R'
'biodyn-hcr.R'
'biodyn-kobe.R'
'biodyn-mng.R'
'biodyn-mse.R'
'biodyn-msy.R'
'biodyn-pella.R'
'biodyn-plots.R'
'biodyn-profile.R'
'biodyns.R'
'combine.R'
'biodyn-cpp.R'\cr
License: \tab GPL-2\cr
LazyLoad: \tab Yes\cr
LazyData: \tab No\cr
Built: \tab R 2.15.3; ; 2013-03-25 16:42:17 UTC; unix\cr
}

Index:
\preformatted{
biodyn                  Biomass Dynamic Model Class
biodyn-class            Class '"biodyn"'
biodyn-package          Biomass Dynamic Model 
biodyns                 Creates a list of biodyn
bmsy                    Biomass at Maximum Sustainable Yield.
calcK                   Calculate Carrying Capacity
computeSP               Calculates surplus production
control<-               control
fit                     fit
fmsy                    Fishing Mortality at Maximum Sustainable Yield.
fwd                     Forward Projections
hcr                     HCR
is.biodyn               Checks class type
msy                     Maximum Sustainable Yield.
pella                   pella
plot                    Create a 'ggplot' plot
plotSP                  Create a 'ggplot' plot
power                   power
setControl<-            control
setParams<-             setParams<-
}

}
\author{
Laurence Kell

Maintainer: Laurence Kell <laurie.kell@iccat.int>
}
\references{

   Fletcher, R.I. 1978. On the restructuring of the Pella-Tomlinson system.
  Fish. Bull. 76:515--521.

   Fox, W.W., Jr. 1970. An exponential surplus-yield model for optimizing
  exploited fish populations. Trans. Am. Fish. Soc. 99:80--88.

   Garrod, D.J. 1969. Empirical assessments of catch effort relationships
  in the North Atlantic cod stock. ICNAF Res. Bull. 6:26--34.

   Gilpin, M.E. and F.J. Ayala. 1973. Global models of growth and
  competition. Proc. Natl. Acad. Sci. USA 70:3590--3593.

   Gompertz, B. 1825. On the nature of the function expressive of the law
  of human mortality, and on a new mode of determining the value of life
  contingencies. Phil. Trans. R. Soc. Lond. 115:513--583.

   Hilborn, R. and C.J. Walters. 1992. Quantitative fisheries stock
  assessment: Choice, dynamics and uncertainty. New York: Chapman and Hall.

  Kingsland, S. 1982. The refractory model: The logistic curve and the
  history of population ecology. Q. Rev. Biol. 57:29--52.

   Laloe, F. 1995. Should surplus production models be fishery description
  tools rather than biological models? Aquat. Living Resour. 8:1--16.

   Pella, J.J. and P.K. Tomlinson. 1969. A generalized stock production
  model. IATTC Bull. 13:421--496.

   Polacheck, T., R. Hilborn, and A.E. Punt. 1993. Fitting surplus
  production models: Comparing methods and measuring uncertainty. Can. J. Fish.
  Aquat. Sci. 50:2597--2607.

   Prager, M.H. 2002. Comparison of logistic and generalized
  surplus-production models applied to swordfish, Xiphias gladius, in the north
  Atlantic Ocean. Fish. Res. 58:41--57.
  
   Quinn, T.J., II and R.B. Deriso. 1999. Quantitative fish dynamics. New
  York: Oxford University Press.

  Schaefer, M.B. 1954. Some aspects of the dynamics of populations
  important to the management of the commercial marine fisheries. IATTC Bull.
  1:27--56.

}

\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
}
