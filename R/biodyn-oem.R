## Observation Error Model

#' sim
#'
#' @description Creates a biodyn object with known properties
#' 
#' @param model character corresponding to model
#' @param params surplus production parameters
#' @param harvest \code{FLQuant} with harvest rate
#' 
#' @param bounds on \code{control}
#' @param ... other arguments
#' 
#' @export
#' @rdname sim
#' 
#' @return biodyn object with simulated time series
#' 
#' @export
#' @examples
#' \dontrun{
#'  bd=sim() 
#'  }
setGeneric('sim',   function(stock,brp,...)     standardGeneric('sim'))
setMethod( 'sim',   signature(stock='missing',brp='missing'),
           function(params=FLPar(r=0.5, k=1000, p=1, b0=1.0),
                    harvest=FLQuant(FLQuant(c(seq(0,1.5,length.out=30), 
                                              rev(seq(0.5,1.5,length.out=15))[-1],
                                              rep(0.5,5)))*biodyn:::fmsy(biodyn(params=params))),
                    bounds =c(0.1,10), ...) {

  args <- list(...)
  
  nyr <- dims(harvest)$year
  object = biodyn(model ='pellat',
                  stock =FLQuant(rep(params['k'], nyr), dimnames=dimnames(harvest)),
                  params=params)
  
  object@control['r',     'val']=params['r']
  object@control['k',     'val']=params['k']
  object@control['p',     'val']=params['p']
  object@control['b0',    'val']=params['b0']
  
  object@control[,'min']=object@control[,'val']*bounds[1]
  object@control[,'max']=object@control[,'val']*bounds[2]
  
  object@control['p', 'phase']=-1
  object@control['b0','phase']=-1
  object@priors[,1]=-1
  
  # Load given slots
  for(i in names(args))
    slot(object, i) <- args[[i]]
  
  object <- fwd(object, harvest=harvest)
  
  return(object)}) 

setMethod( 'sim', signature(stock='FLStock',brp='FLBRP'),function(stock,brp) {
  
  bd=biodyn::biodyn(stock)
  
  params(bd)[dimnames(ctrl)$param]=ctrl[dimnames(ctrl)$param,'val']
  
  bd@priors=prrs
  setParams( bd)=cpue
  setControl(bd)=params(bd)
  bd@control[dimnames(ctrl)$params,'phase'][]=ctrl[dimnames(ctrl)$params,'phase']
  bd@control['q1','phase']=phaseQ
  bd@control['q1','val']  =1
  
  nyr <- dims(harvest)$year
  object = biodyn(model ='pellat',
                  stock =FLQuant(rep(params['k'], nyr), dimnames=dimnames(harvest)),
                  params=params)
  
  object@control['r',     'val']=params['r']
  object@control['k',     'val']=params['k']
  object@control['p',     'val']=params['p']
  object@control['b0',    'val']=params['b0']
  
  object@control[,'min']=object@control[,'val']*bounds[1]
  object@control[,'max']=object@control[,'val']*bounds[2]
  
  object@control['p', 'phase']=-1
  object@control['b0','phase']=-1
  object@priors[,1]=-1
  
  # Load given slots
  for(i in names(args))
    slot(object, i) <- args[[i]]
  
  object <- fwd(object, harvest=harvest)
  
  return(object)})

setGeneric('oem',   function(stock,...)     standardGeneric('oem'))
setMethod( 'oem',   signature(stock='FLStock'),
           function(stock,cv=0.3,fishDepend=FALSE){
  
  nits=max(dims(stock(stock))$iter,dims(catch(stock))$iter)
  rnd=rlnorm(nits,FLQuant(0,dimnames=list(year=dims(stock)$minyear:dims(stock)$maxyear)),cv)
  
  if (fishDepend) 
    cpue=rnd*catch(stock)/fbar(stock)
  else 
    cpue=rnd*computeStock(stock)
  
  cpue})
