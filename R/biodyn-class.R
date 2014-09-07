#' biodyn class
#'
#' @description A class to represent a biomass dynamic stock assessment model. 
#' There are slots for data, parameters and residuals and methods for calculating 
#' reference points and other derived quantities.
#' 
#' @return biodyn object
#' @export
#' @aliases biodyn,character,FLPar-method biodyn,character,missing-method biodyn,factor,FLPar-method
#' biodyn,factor,missing-method biodyn,missing,missing-method
#' 
#' @importFrom plyr ddply ldply
#' @import FLCore
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
    profile       ='data.frame'
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
#' @param object biodyn object or a list of biodyn objects
#' @param ... any other parameter
#' 
#' @return biodyns object
#' @export
#' @rdname biodyns
#' 
#' @aliases biodyns,biodyn-method  biodyns,list-method  biodyns,missing-method
#' 
#' @examples 
#' \dontrun{
#' biodyns(list("a"=bd,"b"=bd))
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