#' biodyn class
#'
#' @description A class to represent a biomass dynamic stock assessment model. 
#' There are slots for data, parameters and residuals and methods for calculating 
#' reference points and other derived quantities.
#' 
#' @return biodyn object
#' @export
#' @aliases biodyn,character,FLPar biodyn,character,missing-method biodyn,factor,FLPar-method
#' biodyn,factor,missing-method biodyn,missing,missing-method
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