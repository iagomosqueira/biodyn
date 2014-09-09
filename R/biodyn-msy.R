utils::globalVariables(c('refJacobian'))
utils::globalVariables(c('jacobian', 'swon'))
utils::globalVariables(c('data.x','data.y','X..x'))

#' msy
#'
#' Calculates MSY given the model parameters, can be done for a biodyn class, or by specifying the model and parameters
#'
#' @param object an object of class \code{biodyn} or a string or factor that species the model
#' @param params an \code{FLPar} object with model parameters
#' @param ... any other parameters
#' 
#' @return an \code{FLPar} object with value(s) of MSY
#' 
#' @seealso \code{\link{bmsy}}, \code{\link{fmsy}} 
#'
#' @aliases msy,biodyn,missing-method  msy,biodyn,numeric-method msy,character,FLPar-method msy,factor,FLPar-method
#'  
#' @export
#' @rdname msy
#'
#' @examples \dontrun{ msy('logistic',FLPar(msy=100,k=500))}
#'   
setGeneric('msy', function(object,params,...) standardGeneric('msy')) 
setMethod('msy',   signature(object='character', params='FLPar'),   function(object=factor(object), params=params)            msyFn(   object, params))
setMethod('msy',   signature(object='factor',    params='FLPar'),   function(object=       object,  params=params)            msyFn(   object, params))
setMethod('msy',   signature(object='biodyn',    params='missing'), function(object,                params=paramas, probs=0.5) msyFn(model(object), params(object),probs=probs))
setMethod('msy',   signature(object='biodyn',    params='numeric'), function(object,                params=params)            msyFn(object, params(object), probs=params))


#' fmsy
#'
#' Calculates $F_{MSY}$ given the model parameters, can be done for a biodyn class, or by specifying the model and parameters
#'
#' @param object an object of class \code{biodyn} or a string or factor that species the model
#' @param params an \code{FLPar} object with model parameters
#' @param ... any other parameters
#'
#' @aliases fmsy,fmsy-method fmsy,biodyn,missing-method  fmsy,biodyn,numeric-method fmsy,character,FLPar-method fmsy,factor,FLPar-method
#' 
#' @return an \code{FLPar} object with value(s) of $F_{MSY}$
#' 
#' @seealso \code{\link{msy}}, \code{\link{bmsy}} 
#' 
#' @export
#' @rdname fmsy
#'
#' @examples
#' \dontrun{
#' fmsy('logistic',FLPar(msy=100,k=500))
#' } 
setGeneric('fmsy', function(object,params,...) standardGeneric('fmsy')) 
setMethod('fmsy',  signature(object='biodyn',    params='missing'), function(object,probs=0.5)                              fmsyFn(model(object), params(object),probs=probs))
setMethod('fmsy',  signature(object='character', params='FLPar'),   function(object=factor(object), params=params)         fmsyFn(object,params))
setMethod('fmsy',  signature(object='factor',    params='FLPar'),   function(object=       object,  params=params)         fmsyFn(object,params))
setMethod('fmsy',  signature(object='biodyn',    params='numeric'), function(object,                params=params)         fmsyFn(object,params(object)))

#' bmsy
#'
#' Calculates $B_{MSY}$ given the model parameters, can be done for a biodyn class, or by specifying the model and parameters
#'
#' @param object an object of class \code{biodyn} or a string or factor that species the model
#' @param params an \code{FLPar} object with model parameters
#' @param ... any other parameters
#'  
#' @aliases bmsy,biodyn,missing-method  bmsy,biodyn,numeric-method bmsy,character,FLPar-method bmsy,factor,FLPar-method
#' 
#' @return an \code{FLPar} object with value(s) of $MSY$, $B_{MSY}$ and $B_{MSY}$
#' 
#' @seealso \code{\link{msy}}, \code{\link{fmsy}} and  \code{\link{bmsy}}
#' 
#' @export
#' @rdname bmsy
#'
#' @examples
#' \dontrun{ 
#' refpts('logistic',FLPar(msy=100,k=500))
#' } 
setGeneric('bmsy', function(object,params,...) standardGeneric('bmsy')) 
setMethod('bmsy',  signature(object='character', params='FLPar'),   function(object=factor(object), params=params)         bmsyFn(  object, params))
setMethod('bmsy',  signature(object='factor',    params='FLPar'),   function(object=       object,  params=params)         bmsyFn(  object, params))
setMethod('bmsy',  signature(object='biodyn',    params='missing'), function(object,probs=0.5)                              bmsyFn(model(object), params(object), probs=probs))
setMethod('bmsy',  signature(object='biodyn',    params='numeric'), function(object,params=params)                         bmsyFn(object, params(object), probs=params))

##############################################################
#' Maximum Sustainable Yield reference points
#'
#' Calculates $MSY$, $B_{MSY}$ and $F_{MSY}$ given the model parameters, can be done for a biodyn class, or by specifying the model and parameters
#'
#' @param object an object of class \code{biodyn} or a string or factor that species the model
#' @param params an \code{FLPar} object with model parameters
#' @param ... any other arguments
#'
#' @return an \code{FLPar} object with value(s) of $F_{MSY}$
#' 
#' @seealso \code{\link{msy}}, \code{\link{bmsy}} 
#' 
#' @export
#' @rdname refpts
#'
#' @aliases refpts refpts-method  refpts,biodyn,missing-method  refpts,factor,FLPar-method   refpts,character,FLPar-method refpts,ANY
#'
#' @examples
#' \dontrun{
#' refpts('logistic',FLPar(msy=100,k=500))
#' } 

setGeneric('refpts',   function(object,params,...) standardGeneric('refpts'))
setMethod('refpts', signature(object='character', params='FLPar'),   function(object=factor(object), params=params)          refptsFn(object, params))
setMethod('refpts', signature(object='factor',    params='FLPar'),   function(object=       object,  params=params)          refptsFn(object, params))
setMethod('refpts', signature(object='biodyn',    params='missing'), function(object)  {  model=model(object)
                                                                                          par  =params(object)
                                                                                          refptsFn(model,par)})
setGeneric('refptSE',  function(object,params,...) standardGeneric('refptSE'))
setMethod('refptSE', signature(object='character', params='FLPar'),   function(object=factor(object), params=params)         refptsFn(object, params))
setMethod('refptSE', signature(object='factor',    params='FLPar'),   function(object=       object,  params=params)         refptsFn(object, params))
setMethod('refptSE', signature(object='biodyn',    params='missing'), function(object)  {  model=model(object)
                                                                                           par  =params(object)
                                                                                           refptsFn(model,par)})


# Fox
msyFox  <- function(params)
  params['r']*(params['k']*exp(-1))*(1-(log(params['k'])-1)/log(params['k']))

# Schaefer
msySchaefer <- function(params)
  params['r']*params['k']/4
msyLogistic <- function(params){
  r=4*params['msy']/params['k']
  r*params['k']/4}

# PellaT
msyPellaT <- function(params)
  params['r']*params['k']*(1/(1+params['p']))^(1/params['p']+1)
msyGenfit <- function(params)
  params['r']*params['k']*(1/(1+params['p']))^(1/params['p']+1)

# Shepherd
msyShepherd<-function(params) {
  aPrime<-params['r']/params['m'] - 1
  Bmax  <-params['k']*aPrime
  .bmsy <- 0 #bmsy('shepherd',param)
  
  aPrime*params['m']*.bmsy*(1-.bmsy/Bmax)/(1+aPrime)^.5}

# Gulland
msyGulland  <- function(params)
  (params['r']*params['k']^2)/4

msyFletcher <- function(params)
  params['msy']

# Fox
bmsyFox  <- function(params)
  params['k']*exp(-1)
# Schaefer
bmsySchaefer <- function(params)
  params['k']/2
bmsyLogistic <- function(params)
  params['k']/2
# PellaT
bmsyPellaT <- function(params)
  params['k']*(1/(1+params['p']))^(1/params['p'])
bmsyGenfit <- function(params)
  params['k']*(1/(1+params['p']))^(1/params['p'])
# Shepherd
bmsyShepherd <- function(params) {
  aPrime <- params['r']/params['m'] - 1
  Bmax  <- params['k']*aPrime
  
  Bmax*((1+aPrime)^.5-1)/aPrime}

# Gulland
bmsyGulland <-function(params)
  params['k']/2
# Fletcher
bmsyFletcher <- function(params)
  params['k']*(1/(params['p']+1)^(1/(params['p'])))

fmsyFn=function(object,params,probs=0.5){
  
  if (probs!=0.5){
    res=qnorm(probs, object@mng['fmsy','hat'],object@mng['fmsy','sd'])
    return(res)}
  
  object=tolower(object)
  fmsyPellaT  <-function(params) params['r']*(1/(1+params['p']))
  fmsyFox     <-function(params) params['r']*(1-(log(params['k'])-1)/log(params['k']))
  fmsySchaefer<-function(params) params['r']/2
  fmsyShepherd<-function(params) msyShepherd(params)/bmsyShepherd(params)
  fmsyGulland <-function(params) params['r']*params['k']/2
  fmsyFletcher<-function(params) msyFletcher(params)/bmsyFletcher(params)
  fmsyLogistic<-function(params) {r=4*params['msy']/params['k']; r/2}
  
  res<-switch(as.character(object),
              
              fox     =fmsyFox(     params),
              schaefer=fmsySchaefer(params),
              gulland =fmsyGulland( params),
              fletcher=fmsyFletcher(params),
              pellat  =fmsyPellaT(  params),
              shepherd=fmsyShepherd(params),
              logistic=fmsyLogistic(params),
              genfit  =fmsyPellaT(params))
  
  dimnames(res)$params='fmsy'
  
  return(res)}

msyFn=function(object,params,probs=0.5) {
   
  if (probs!=0.5){
    res=qnorm(probs, object@mng['msy','hat'],object@mng['msy','sd'])
    return(res)}
  
  
  object=tolower(object)
  res<-switch(object,
              fox      = msyFox(params),
              schaefer = msySchaefer(params),
              gulland  = msyGulland(params),
              fletcher = msyFletcher(params),
              pellat   = msyPellaT(params),
              shepherd = msyShepherd(params),
              logistic = msyLogistic(params),
              genfit   = msyPellaT(params))
  
  dimnames(res)$params='msy'
  
  return(res)}

bmsyFn=function(object,params,probs=0.5) {
  
  if (probs!=0.5){    
      res=qnorm(probs, object@mng['bmsy','hat'],object@mng['bmsy','sd'])
      return(res)}
  
  object=tolower(object)
  res<-switch(object,
              fox     =bmsyFox(params) ,
              schaefer=bmsySchaefer(params),
              gulland =bmsyGulland(params),
              fletcher=bmsyFletcher(params),
              pellat  =bmsyPellaT(params),
              shepherd=bmsyShepherd(params),
              logistic=bmsySchaefer(params),
              genfit  =bmsyPellaT(params))
  
  dimnames(res)$params='bmsy'
  
  return(res)}

##############################################################
#' Calculate Carrying Capacity
#'
#' Calculates $k$ given msy, r and K for a Pella-Tomlinson biomass dynamic model
#'
#' @param msy a guess for $MSY$
#' @param r a guess for $r$ the population growth rate
#' @param p a guess for $p$ the shape parameter
#' @param params provide $r$ and $p$ as \code{FLPar}
#'
#' @return an \code{FLPar} object with an estimate for $k$
#' 
#' @seealso \code{\link{msy}} and \code{\link{bmsy}} 
#' 
#' @export
#' @rdname K
#'
#' @examples
#' \dontrun{
#'  K(5000,r=.6,p=1.0)
#'  }
K <- function(msy,r=.6,p=1,params=FLPar(r=r,p=p)){
  res=msy/(params['r']*(1/(1+params['p']))^(1/params['p']+1))
  
  dimnames(res)$params='k'
  
  return(res)} 

refptsFn=function(model,params){
  model=tolower(model)
  dmns<-dimnames(params)
  names(dmns)[1]<-'refpts'
  dmns[[1]]<-c('msy','fmsy','bmsy')
  res<-FLPar(as.numeric(NA),dimnames=dmns)
  
  obj=as.character(model)
  
  res[ 'msy']<- msyFn(model, params)
  res['fmsy']<-fmsyFn(model, params)
  res['bmsy']<-bmsyFn(model, params)
  
  return(res)}

fnB=function(x,bd,year=range(bd)['maxyear']){
  if (is.numeric(year)) year=ac(year)
  params(bd)[c('r','k')]=x
  stock(bd[,year])/refpts(bd)['bmsy']}

fnF=function(x,bd,year=range(bd)['maxyear']){
  if (is.numeric(year)) year=ac(year)
  params(bd)[c('r','k')]=x
  harvest(bd[,'2008'])/refpts(bd)['fmsy']}

fnY=function(x,bd,year=range(bd)['maxyear']){
  if (is.numeric(year)) year=ac(year)
  params(bd)[c('r','k')]=x
  catch(bd[,'2008'])/refpts(bd)['msy']}

fnJ=function(x,bd,year=range(bd)['maxyear']){
  if (is.numeric(year)) year=ac(year)
  bd@params[activeParams(bd)]=x
  c(stock  =stock(  bd[,'2008'])/refpts(bd)['bmsy'],
    harvest=harvest(bd[,'2008'])/refpts(bd)['fmsy'],
    catch  =catch(  bd[,'2008'])/refpts(bd)[ 'msy'])}  


refJacobian=function(bd,year=range(bd)['maxyear']){

  x0=bd@params[activeParams(bd)]
  
  res <- jacobian(func=fnJ, x=x0, bd=swon)
  
  nms=dimnames(bd@params[activeParams(bd)])[[1]]
  
  res=FLPar(array(t(res),dim=c(length(nms),3,1),dimnames=list('params'=nms,' '=c('stock','harvest','catch'),iter=1)))
  
  return(res)}

relVar=function(bd,year=range(bd)['maxyear']){
  
  nms=dimnames(bd@params[activeParams(bd)])[[1]]
  res=refJacobian(bd,year=year)
  cov=bd@vcov[nms,nms]
  
  
  res.=as.data.frame(res)
  
  FLCore::iter(cov,1)[lower.tri(iter(cov,1)[drop=T])] =NA
  cov.=as.data.frame(cov)
  cov.=cov.[!is.na(cov.$data),]
  
  vars=merge(cov., merge(res.,res.,by='params',all=T),by='params',all=T)
  
  vars=ddply(transform(vars,var=data.x*data.y*data)[,c('X..x','iter','var')],.(X..x,iter), function(x) data.frame(var=sum(x$var)))
  
  return(res)}  



