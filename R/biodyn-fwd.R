# fwd.R - 
# biodyn/R/fwd.R

# Copyontright 2003-2007 FLR Team. Distributed under the GPL 2 or later
# Maintainer: Iago Mosqueira, JRC
# $Id:  $

utils::globalVariables('finite')

prdFn=function(mdl,params,biomass=0) {
  if (!is.FLQuant(biomass)) biomass=FLQuant(biomass)  
  
  if (dims(params)$iter==1 & dims(biomass)$iter >1) params =propagate(params,dims(biomass)$iter)
  if (dims(params)$iter> 1 & dims(biomass)$iter==1) biomass=propagate(biomass,dims(params)$iter)
  
  mdl=tolower(as.character(mdl))
  
  fox <-function(biomass, params)
    params['r']*biomass*(1-log(biomass)/log(params['k']))
  
  schaefer <- function(biomass, params)
    params['r']*biomass*(1-biomass/params['k'])
  
  pellat <- function(biomass, params){
    a=sweep(biomass,6,params['r']/params['p'],'*')
    b=sweep(biomass,6,params['k'],'/')
    c=sweep(b,      6,params['p'],'^')
    a*(1-c)}
  
  shepherd <- function(biomass,params)
    params['r']*biomass/(1+biomass/params['k'])-params['m']*biomass
  
  gulland <- function(biomass,params)
    params['r']*biomass*(params['k']-biomass)
  
  fletcher <- function(biomass,params) {
    params['p']=params['p']+1
    lambda <- (params['p']^(params['p']/(params['p']-1)))/(params['p']-1)
    lambda*params['msy']*(biomass/params['k'])-lambda*params['msy']*(biomass/params['k'])^params['p']
  }
  
  logistic <- function(biomass, params){
    r=4*params['msy']/params['k']
    r*biomass%*%(1-biomass%/%params['k'])}
  
  genfit <- function(biomass, params)
    params['r']*biomass*(1-biomass/params['k'])
  
  res <- switch(mdl,
                fox     =fox(     biomass,params),
                schaefer=schaefer(biomass,params),
                gulland =gulland( biomass,params),
                fletcher=fletcher(biomass,params),
                pellat  =pellat(  biomass,params=params),
                shepherd=shepherd(biomass,params),
                genfit  =pellat(  biomass,params),
                logistic=logistic(biomass,params))
  
  return(res)}

# calcs a relative value
relFn=function(val,ratio,lag=0){

  if (lag<=0) return(ratio)
  
  yrs=ac(as.numeric(dimnames(ratio)$year)-lag)
  yrs=yrs[yrs %in% dimnames(val)$year]
  yr2=ac(as.numeric(yrs)+lag)
  res=ratio[,yr2]*val[,yrs]
  
  return(res)}

## sets a minimum F
minFn=function(val,bnd){
  val[is.na(val)]=0
  
  res=val
  res[val<bnd]=bnd[val<bnd]
  
  return(res)}

## sets a maximum F
maxFn=function(val,bnd){
  val[is.na(val)]=0
  
  res=val
  res[val>bnd]=bnd[val>bnd]
 
  return(res)}

## calculates inter-annual bounds
iavFn=function(val,bnd,lag=1){
  if (all(!is.finite(bnd)) | lag<1) stop('cant have lags < 1')
  
  ref=apply(val[,dims(val)[2]-lag],2:6,mean)
  
  bnd=val[,dims(val)[2]]*(1-bnd[1])    
  if (finite(bnd[1])) val[,dims(val)[2]][ref<bnd[1]]=ref[ref<bnd[1]]
  
  bnd=val[,dims(val)[2]]*(1+bnd[2])    
  if (finite(bnd[2])) val[,dims(val)[2]][ref>bnd[2]]=ref[ref>bnd[2]]
 
  return(val[,dims(val)[2]])}


#val=FLQuant(rlnorm(200),dimnames=list(age=1:5,year=1:10,iter=1:4))
#bnd=FLQuant(1,          dimnames=list(age=1:5,year=1:10,iter=1:4))


#' fwd
#'
#' @description Projects a \code{biodyn} object (i.e. a biomass dynamic model)
#' for a given future catch, harvest rate, or stock biomass. Only
#' one of these has to be supplied as an argument.
#'
#' @param object an object of class \code{biodyn} or  \code{biodyns}
#' @param ctrl missing
# @param catch   an \code{FLQuant} or \code{FLQuants} containing future catches
# @param harvest an \code{FLQuant} or \code{FLQuants} containing future harvest
# @param stock   an \code{FLQuant} or \code{FLQuants} containing future stock    
# @param minF    minimum harvest, \code{FLQuant} or single numeric value    
# @param maxF    maximum harvest, \code{FLQuant} or single numeric value
# @param pe      process error term, an \code{FLQuant}
# @param peMult  logical, i.e. by default equals true so multiplicative, otherwise additive
#' @param ... any other parameters
#'
#' @aliases fwd-method fwd,biodyn,FLQuants-method  fwd,biodyn,missing-method 
#' 
#' @export
#' @rdname fwd
#'
#' @examples
#' \dontrun{
#'  bd=biodyn('logistic',FLPar(msy=1000,k=50000,b0=1),stock=FLQuant(50000,dimnames=list(year=1:50)))
#'  bd=fwd(bd,catch=FLQuant(500,dimnames=list(year=1:50)))
#'  }
setGeneric('fwd', function(object,ctrl,...) standardGeneric('fwd'))
setMethod( 'fwd', signature(object='biodyn',ctrl='missing'),
   function(object, catch  =NULL, 
                    harvest=NULL, 
                    stock  =NULL, 
                    hcr    =NULL, 
                    pe     =NULL, peMult=TRUE,
                    minF   =0,    maxF  =2,
                    # these next arguments are not currently used, as 
                    # the code havnt been enabled
                    bounds =list(catch=c(Inf,Inf)),
                    lag    =0,
                    end    =NULL,
                    ...) {
     
     ## target arg is an FLQuant
     
     if ('FLQuant' %in% class(stock)   |
         'FLQuant' %in% class(harvest) |
         'FLQuant' %in% class(catch))
        res=fwdFn(object,ctrl=ctrl,
              catch,harvest,stock,pe,peMult,minF,maxF,bounds,lag,end,...) 
     else if ('FLQuants' %in% class(stock))  
        res=biodyns(llply(stock, function(x) fwdFn(object,ctrl=ctrl,
                 catch,harvest,x,pe,peMult,minF,maxF,bounds,lag,end,...)))
      else if ('FLQuants' %in% class(harvest))  
        res=biodyns(llply(harvest, function(x) fwdFn(object,ctrl=ctrl,
                 catch,x,stock,pe,peMult,minF,maxF,bounds,lag,end,...))) 
      else if ('FLQuants' %in% class(catch))  
        res=biodyns(llply(catch, function(x) fwdFn(object,ctrl=ctrl,
                x,harvest,stock,pe,peMult,minF,maxF,bounds,lag,end,...)))
                        
     res})
     
fwdFn=function(object, 
               catch, 
               harvest, 
               stock, 
               pe, peMult,
               minF,    maxF,
               ...){
        
      lag=0
      object@stock=FLQuant(object@stock,quant=names(object@catch)[1])   
    
      ## catch, harvest or stock?
      ctcTrgt=FALSE
      hvtTrgt=FALSE
      stkTrgt=FALSE
      hcrTrgt=FALSE
      
      if (!is.null(catch))   ctcTrgt=TRUE 
      if (!is.null(harvest)) hvtTrgt=TRUE
      if (!is.null(stock))   stkTrgt=TRUE
      ## hcr option not implemented, see hcr method
      
      if (hvtTrgt) f=harvest
      
      if(!ctcTrgt & hvtTrgt & stkTrgt & hcrTrgt)
          stop('must supply catch, harvest or stock as a target')
    
      if (ctcTrgt) if (dims(object@stock)$maxyear   < dims(  catch)$maxyear) object=window(object,end=dims(  catch)$maxyear)
      if (hvtTrgt) if (dims(object@stock)$maxyear-1 < dims(harvest)$maxyear) object=window(object,end=dims(harvest)$maxyear+1)
      if (stkTrgt) if (dims(object@stock)$maxyear   < dims(  stock)$maxyear) object=window(object,end=dims(  stock)$maxyear)
       
      if (stkTrgt) catch=stock*0
     
      ## check year range
      if (ctcTrgt | stkTrgt) {
        if (!(all(dimnames(catch)$year %in% dimnames(object@catch)$year)))
             object = window(object,end=dims(catch)$maxyear)
          if (dims(object@catch)$iter==1 & dims(catch)$iter>1) object@catch=propagate(object@catch, dims(catch)$iter)
          object@catch[,dimnames(catch)$year] <- catch
          yrs <- dimnames(catch)$year
      } else if (hvtTrgt) {
          if (!(all(dimnames(harvest)$year %in% dimnames(object@stock)$year))){
            
            stop('years in harvest & stock dont match')}
          yrs <- dimnames(harvest)$year
      } else if (hcrTrgt) {
         yrs=ac(dims(object@stock)$maxyear:end) 
         object=window(object,end=end)
      }  
     
      if (stkTrgt) yrs = yrs[-length(yrs)]
      ## B0 in year 1?
      if (as.numeric(yrs[1]) == range(object,'minyear')){
         if (!('year' %in% names(dimnames(params(object)['k']))))  
           object@stock[,ac(range(object,'minyear'))] = params(object)['k'] * params(object)['b0'] else
           object@stock[,ac(range(object,'minyear'))] = params(object)['k',ac(range(object,'minyear'))] * params(object)['b0',ac(range(object,'minyear'))]          
         }
      
      ## maxyear
      if (max(as.numeric(yrs)) == range(object,'maxyear'))
         object@stock <- window(object@stock,end=range(object,'maxyear')+1)
    
      ## niters
      nits=dims(object)$iter
      if (hvtTrgt) nits=max(nits,dims(harvest)$iter)
      if (ctcTrgt) nits=max(nits,dims(  catch)$iter)
      if (stkTrgt) nits=max(nits,dims(  stock)$iter)
      if (!is.null(pe)) nits=max(nits,dims(pe)$iter)
      
      if (hvtTrgt) nits=max(nits,dims(harvest)$iter) else
      if (ctcTrgt) nits=max(nits,dims(catch  )$iter) else
      if (stkTrgt) nits=max(nits,dims(stock  )$iter) 
      if (nits>1){ 
         if(dim(object@catch)[6]==1) object@catch =propagate(object@catch,nits)
         if(dim(object@stock)[6]==1) object@stock =propagate(object@stock,nits)
         if (hvtTrgt) if(dim(harvest)[6]==1) harvest=propagate(harvest,nits)
       
         if (dims(params(object))$iter==1) params(object)=propagate(params(object),nits)
         #if (!is.null(pe))                 pe            =propagate(pe            ,nits)
         } 
      
      ## projections
      if (!hvtTrgt) harvest=harvest(object)
      for(y in as.numeric(yrs)) {
    
         ## sp & process error
         if (!is.null(pe)) {    
            if (peMult) sp.=computePrd(object,object@stock[, ac(y)])%*%pe[, ac(y)] 
            else        sp.=computePrd(object,object@stock[, ac(y)])%+%pe[, ac(y)]
         } else sp.=computePrd(object,object@stock[, ac(y)])
         #} else sp.=computePrd(object,object@stock[, ac(y)]*(1-ptYr)+object@stock[, ac(y+1)]*(ptYr))
         
          ## targets, if lag<0 then the targets are relative 
          if (hcrTrgt)
              harvest[,ac(y)]=hcr(window(object,end=y))
          if (hvtTrgt | hcrTrgt){
              object@catch[,ac(y)]=object@stock[,ac(y)]*relFn(harvest(object),harvest[,ac(y)],lag)
          } else if (ctcTrgt){           
              object@catch[,ac(y)]=relFn(object@catch,catch[,ac(y)],lag)
          } else {
              object@catch[,ac(y)]=relFn(object@stock,stock[,ac(y+1)],lag) - object@stock[,ac(y)] - sp.}
         
         object@stock[,ac(y+1)] = object@stock[,ac(y)] - object@catch[,ac(y)] + sp.
    
    ## Not implemented     
    ### bounds
    #      if ('catch' %in% bounds){ 
    #           object@catch[,y]=iavFn(window(object@catch,end=y),bnd,lag=lag)
    #           object@stock[,ac(y+1)] = object@stock[,ac(y)] - object@catch[,ac(y)] + sp.
    #           }
        
          ### min/max
          harvest[,ac(y)]       =maxFn(minFn(harvest(object)[,ac(y)],minF),maxF)
          object@catch[,ac(y)]  =object@stock[,ac(y)]*harvest[,ac(y)]
          object@stock[,ac(y+1)]=object@stock[,ac(y)] - object@catch[,ac(y)] + sp.
          }
      
        object@stock[stock(object) < 0] = 0.001
        object@catch[catch(object) < 0] = 0.0000001
    
        return(object)}


setMethod('fwd', signature(object='biodyn',ctrl='FLQuants'),
  function(object, ctrl, pe=NULL, peMult=TRUE,minF=0,maxF=2,lag=0,
           bounds=list(catch=c(Inf,Inf)),...) {
    
  res=mlply(seq(length(names(ctrl))),
      function(x,object,ctrl,pe,peMult,minF,maxF,lag,bounds){
        if (names(ctrl)[x]=='catch')  return(fwd(object,catch  =ctrl[[x]],pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds))
        if (names(ctrl)[x]=='harvest')return(fwd(object,harvest=ctrl[[x]],pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds))
        if (names(ctrl)[x]=='stock')  return(fwd(object,stock  =ctrl[[x]],pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds))},
       object=object,ctrl=ctrl,pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds)
  
  names(res)=names(ctrl)
        
  return(biodyns(res))})

## this is wrapper for a HCR projection, but now being done by the HCR method
# setMethod('fwd', signature(object='biodyn',ctrl='list'),
#   function(object, ctrl, pe=NULL, peMult=TRUE,minF=0,maxF=2,lag=0,
#            bounds=list(catch=c(Inf,Inf)),end=range(object,'maxyear')+15,...) {
#     
#   res=mlply(seq(length(names(ctrl))),
#       function(x,object,ctrl,pe,peMult,minF,maxF,lag,bounds,end){
#         return(fwd(object,hcr=ctrl[[x]],pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds,end=end))},
#        object=object,ctrl=ctrl,pe=pe,peMult=peMult,minF=minF,maxF=maxF,lag=lag,bounds=bounds,end=end)
#   
#   names(res)=names(ctrl)
#         
#   return(biodyns(res))})
      

    
