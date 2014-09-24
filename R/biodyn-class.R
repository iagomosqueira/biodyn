#' biodyn class
#'
#' @description A class that implement a biomass dynamic stock assessment model. 
#' 
#' @details 
#' The Class is intended to be used as part of an MSE and includes methods for diagnostics, calculating  reference points and other quantities used when providing management advice.
#' 
#' @slot name    {A \code{character} with the name of the stock}
#' @slot desc    {A \code{character} providing a fuller description of the object}       
#' @slot range   {A \code{numeric} vector containing the quant and year ranges}
#' @slot model   {A \code{factor} giving name of production function, for now this is only `pellat`}
#' @slot catch   {An \code{FLQuant} with total catch by year}        
#' @slot stock   {An \code{FLQuant} which will hold the estimated stock by yaer}       
#' @slot control {An \code{FLPar} which sets initial guess (val) and bounds (min and max) for each parameter. The phase allows a parameter to be fixed if less <0 and for paramters to be estimated sequentially}       
#' @slot priors  {An \code{array} which sets penalties for parameters}         
#' @slot params  {An \code{FLPar} with parameter estmates}
#' @slot vcov    {An \code{FLPar} with the covariance matrix of the parameter estimates} 
#' @slot hessian {An \code{FLPar} with the hessian of the estimated parameters} 
#' @slot mng     {\code{FLPar} with derived quatities of management interest}
#' @slot mngVcov {An \code{FLPar} with the variance matrix of management quanties}
#' @slot diags   {A \code{data.frame} with residuals and covariates from fit of CPUE to stock }     
#' @slot objFn   {\code{FLPar} with objective function} 
#' @slot ll      {\code{FLPar} with negative log likelihood by data component}
#' @slot profile {\code{data.frame} not yet implemented} 
#' @slot hst     {\code{data.frame} not yet implemented} 
#' 
#' All slots in the class have accessor and replacement methods that provide validation and protection of their data.
#' 
#' @export
#' 
#' @importFrom plyr ddply ldply laply mdply maply alply .
#' @import reshape 
#' @importFrom stringr str_trim
#' @import FLCore
#' @import methods
#' 
#' @examples
#' \dontrun{biodyn()}
 setClass('biodyn', representation(
    'FLComp',
    model         ='factor',
    catch         ='FLQuant',
    stock         ='FLQuant',
    diags         ='data.frame',
    params        ='FLPar',
    control       ='FLPar',
    priors        ='array',
    vcov          ='FLPar',
    hessian       ='FLPar',
    objFn         ='FLPar',
    ll            ='FLPar',
    mng           ='FLPar',
    mngVcov       ='FLPar',
    profile       ='data.frame',
    hst           ='data.frame'
    ),
  prototype(
    range       =unlist(list(minyear=as.numeric(NA), maxyear=as.numeric(NA))),
    catch       =FLQuant(),
    stock       =FLQuant(),
    model       =models[3],
    params      =FLPar(c(.5,as.numeric(NA),2,1,as.numeric(NA),as.numeric(NA)),                            dimnames=list(params=c('r','k','p','b0'),iter=1)),
    control     =FLPar(array(rep(c(1,as.numeric(NA),as.numeric(NA),as.numeric(NA)),each=4), dim=c(4,4,1), dimnames=list(params=c('r','k','p','b0'),option=c('phase','min','val','max'),iter=1))),
    priors      =array(rep(c(0,0,0.3,1),       each=7), dim=c(7,4),   dimnames=list(params=c('r','k','p','b0','msy','bmsy','fmsy'),c('weight','a','b','type'))),
    vcov        =FLPar(array(as.numeric(NA), dim=c(4,4,1), dimnames=list(params=c('r','k','p','b0'),params=c('r','k','p','b0'),iter=1))),
    hessian     =FLPar(array(as.numeric(NA), dim=c(4,4,1), dimnames=list(params=c('r','k','p','b0'),params=c('r','k','p','b0'),iter=1))),
    objFn       =FLPar(array(as.numeric(NA),dim=c(2,1),dimnames=list('value'=c('ll','rss'),iter=1)))),
	validity=validity) 

#' biodyns
#' @description Create a list with biodyn objects
#' @name biodyns
#' @param object \code{biodyn} object or a \code{list} of \code{biodyn} objects
#' @param ... additional \code{biodyn} objects
#' 
#' @return \code{biodyns} object
#' @export
#' @rdname biodyns
#' 
#' @aliases biodyns,biodyn-method  biodyns,list-method  biodyns,missing-method
#' 
#' @examples 
#' \dontrun{
#' biodyns(biodyn())
#' }

setGeneric('biodyns', function(object, ...) standardGeneric('biodyns'))
biodyns <- setClass('biodyns', contains='FLComps',
                    validity=function(object) {
                      # All items are biodyn
                      if(!all(unlist(lapply(object, is, 'biodyn'))))
                        return('Components must be biodyn')  
                      
                      return(TRUE)})

setGeneric('biodyns', function(object, ...)
  standardGeneric('biodyns'))

# constructor
setMethod('biodyns', signature(object='biodyn'), function(object, ...) {
  lst <- c(object, list(...))
  biodyns(lst)
})

setMethod('biodyns', signature(object='missing'),
          function(...) {
            # empty
            if(missing(...)){
              new('biodyns')
              # or not
            } else {
              args <- list(...)
              object <- args[!names(args)%in%c('names', 'desc', 'lock')]
              args <- args[!names(args)%in%names(object)]
              do.call('biodyns',  c(list(object=object), args))
            }
          }
)

setMethod('biodyns', signature(object='list'),
          function(object, ...) {
            
            args <- list(...)
            
            # names in args, ... 
            if('names' %in% names(args)) {
              names <- args[['names']]
            } else {
              # ... or in object,
              if(!is.null(names(object))) {
                names <- names(object)
                # ... or in elements, ...
              } else {
                names <- unlist(lapply(object, name))
                # ... or 1:n
                idx <- names == 'NA' | names == ''
                if(any(idx))
                  names[idx] <- as.character(length(names))[idx]
              }
            }
            
            # desc & lock
            args <- c(list(Class='biodyns', .Data=object, names=names),
                      args[!names(args)%in%'names'])
            
            return(
              do.call('new', args)
            )
            
          }) 