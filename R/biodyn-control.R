utils::globalVariables(c('ldply','melt','variable'))

#if (!isGeneric('control'))    
  setGeneric('control',     function(object,method,...) standardGeneric('control'))
#if (!isGeneric('control<-'))  
  setGeneric('control<-',   function(object,value)      standardGeneric('control<-'))

#' setParams<-
#'
#' Add parameters for catchability and sigma in the param slot given an index or indices
#' of stock abundance as \code{FLQuant} or \code{FLQuants}. 
#'
#' @param object an object of class \code{biodyn}
#' @param value \code{params}
#'
#' @export
#' @rdname setParams
#' 
#' @aliases setParams<-,biodyn,FLPar-method  setParams<-,biodyn,FLQuant-method setParams<-,biodyn,FLQuants-method  setParams<-,biodyn,data.frame-method
#' 
#' @examples
#' \dontrun{
#' bd=simBiodyn()
#' cpue=stock(bd)
#' setParams(bd) <-cpue
#' setParams(bd)}
#' 
setGeneric('setParams<-', function(object,value,...)  standardGeneric('setParams<-'))

setMethod('setParams<-', signature(object='biodyn',value='data.frame'), function(object,value) {
  nms=c(modelParams(as.character(object@model)),'b0')
  object@params=object@params[nms]

  object@params =setQ(object,value)
  
  return(object)})

setMethod('setParams<-', signature(object='biodyn',value='FLPar'), function(object,value) {
  object@params=value 
  
  return(object)})

setMethod('setParams<-', signature(object='biodyn',value='FLQuant'), function(object,value) {

  nms=c(modelParams(tolower(as.character(object@model))),'b0')
  #value=FLCore::apply(value,2,mean)
  object@params =setQ(object,value)
  
  return(object)})

setMethod('setParams<-', signature(object='biodyn',value='FLQuants'), function(object,value,msy=TRUE) {
  nms=c(modelParams(as.character(object@model)),'b0')
  object@params=object@params[nms]
    
  object@params =setQ(FLCore::iter(object,1),FLQuants(lapply(value, function(x) FLCore::apply(x,2,mean,na.rm=T))))
  
  nms=dimnames(params(object))$param
  n  = as.numeric(summary(substr(nms,1,1)=='q')['TRUE'])
  nms[substr(nms,1,1)=='q']    =paste('q',    seq(n),sep='')
  nms[substr(nms,1,5)=='sigma']=paste('sigma',seq(n),sep='')
  
  dimnames(params(object))$params=nms
  
  return(object)})

getP<-function(bmsy,k,p=c(0.001,5)){
  optimise(function(p,bmsy,k) 
    (bmsy-k*(1/(1+p))^(1/p))^2, p, bmsy=bmsy,k=k)$minimum}

getRK<-function(msy,bmsy,k,p=c(0.0001,5)){
  p=getP(bmsy,k,p)
  k=bmsy/((1/(1+p))^(1/p))
  r=msy/(k*(1/(1+p))^(1/p+1))
  
  FLPar(c(r=r,k=k,p=p))}

#getRK(100,1000,2000)

setMethod('setParams<-', signature(object='biodyn',value='FLBRP'), function(object,value,msy=TRUE) {
  
  if (msy){
    params(object)[c("r","k","p")]=getRK(value@refpts["msy",   "yield"],
                                         value@refpts["msy",   "biomass"],
                                         value@refpts["virgin","biomass"])
    object@priors[c("r","k","p"),"a"]=params(object)[c("r","k","p")]
  }else{
    require(popbio)
    params(object)["k"]=value@refpts["virgin","biomass"]
    params(object)["p"]=getP(c(value@refpts["msy","biomass"]),
                             c(params(object)["k"]),c(0.001,5))
    params(object)["r"]=log(lambda(leslie(eql,fbar=c(value@refpts["crash","harvest"]))))
    }
      
  object=fwd(object,catch=catch(object)[,-1])
  
  return(object)})


#' setControl<-
#'
#' @description Sets the control slot in a biodyn object given the parameters in the \code{params}
#' slot. The starting values \code{val} are set to those in \code{params} and the \code{min} and
#' \code{max} bounds to .1 and 10 times these.
#'
#' @param object an object of class \code{biodyn}
#' @param value  \code{params} object
#' @param ... any other parameter
#'
#' @seealso \code{\link{controlFn}}
#' 
#' @export

#' @rdname setControl
#'
#' @aliases setControl<-,biodyn,FLPar-method  setControl<-,biodyn,FLQuant-method  setControl<-,biodyn,FLQuants-method
#' 
#' @examples
#' \dontrun{
#' bd=simBiodyn()
#' setControl(bd)=params(bd)
#' }
#'  
setGeneric('setControl<-',  function(object,...,value) standardGeneric('setControl<-'))
setMethod('setControl<-', signature(object='biodyn',value='FLPar'), function(object,value,min=0.1,max=10.0) {
  
  if (dims(value)$iter>1 & dims(object@control)$iter==1)
    object@control=propagate(control(object),dims(value)$iter)

  ctr=object@control
  nms=dimnames(object@params)$params
  
  object@control=FLPar(array(rep(c(1,NA,NA,NA),each=length(nms)), dim=c(length(nms),4,dims(value)$iter), dimnames=list(params=nms,option=c('phase','min','val','max'),iter=seq(dims(value)$iter))))
  nms.=nms[nms %in% dimnames(ctr)$params]
  
  object@control[nms.,'phase']=ctr[nms.,"phase"]
  nms.=nms[nms %in% dimnames(object@params)$params]

  object@control[,"val"]=object@params[nms.]
  object@control[,"min"]=object@params[nms.]*min
  object@control[,"max"]=object@params[nms.]*max
  
  if (!is.na(any(value[nms]<0)) & any(value[nms]<0))
    object@control[nms[value[nms]<0],c('min','max')]=object@control[nms[value[nms]<0],c('max','min')]

  prr=object@priors
  nms=c(nms[!substr(dimnames(object@params)$params,1,1) %in% c('s','q')],
        'msy','bmsy','fmsy',
        nms[substr(dimnames(object@params)$params,1,1) %in% c('s','q')])
  object@priors=array(rep(c(0,0,0.3,1),each=length(nms)), dim=c(length(nms),4),   dimnames=list(params=nms,c('weight','a','b','type')))
  nms=dimnames(prr)$params[dimnames(prr)$params %in%  dimnames(object@priors)$params]
  
  object@priors[nms,]=prr[nms,]
  
  return(object)})

.calcSigma <- function(obs,hat=rep(0,length(obs)),na.rm=T){
  n  =length(obs[!is.na(obs+hat)])
  SS =sum((obs-hat)^2,na.rm=na.rm)

  return((SS/n)^.5)}

calcB0<-function(index,q,k,nyrB0=3,error='log'){
  if (is.null(nyrB0)) return(params['b0'])
  
  if (error=='log'){
    t.<-sweep(log(index[,1:nyrB0,,,,,drop=FALSE]),c(1,6),q,'/')
    return(qmax(qmin(exp(apply(t.,c(1,6),mean))/k,1),0))}
  if (error=='normal'){
    t.<-sweep(index[,1:nyrB0,,,,,drop=FALSE],c(1,6),q,'/')
    return(qmax(qmin(apply(t.,c(1,6),mean)/k,1),0))}       
}

calcQ<-function(stock,index,error='log',na.rm=T){

  stock<-(stock[-length(stock)]+stock[-1])/2
  n    <-length(stock)
  index<-index[seq(n)]
  if (na.rm)
    n=length(seq(n)[!is.na(index+stock)])

  res=switch(error,
             normal={q    =sum(stock*index, na.rm=T)/sum(stock*stock, na.rm=na.rm)
                     sigma=.calcSigma(index/(q*stock))
                     data.frame(q=q,sigma=sigma)
             },
             log   ={q    =exp(sum(log(index)-log(stock), na.rm=na.rm)/n)
                     sigma=.calcSigma(log(index),log(q*stock))
                     data.frame(q=q,sigma=sigma)},
             cv   ={res   <-sum(index/stock)
                    sigma2<-.calcSigma(res,na.rm=na.rm)
                    q     <-(-res+(res^2+4*length(index)*sigma2*sum((index/stock)^2)))/(2*length(index)*sigma2)
                    data.frame(q=q,sigma=sigma)})
  
  return(res)}

setQ=function(object,value,error='log'){
  fn=function(value,stock){
    if (dims(value)$iter==1)
      dimnames(value)$iter=1
    if (dims(stock)$iter==1)
      dimnames(stock)$iter=1
    
    if (dims(stock)$iter==1 & dims(value)$iter>1)
      stock=propagate(stock,dims(value)$iter)
    
    model.frame(mcf(FLQuants(stock=stock,value=value)))}

  res=switch(is(value)[1],
             FLQuant   ={res=fn(value,stock(object));data.frame(name=1,res)},
             #FLQuants  =ldply(value, fn, model.frame(mcf(FLQuants(stock=stock,value=x))),stock=stock(object)),
             FLQuants  =ldply(value, function(x,stock) fn(x,stock), stock=stock(object)),
             data.frame=merge(model.frame(FLQuants('stock'=stock(object))),value,by='year',all=T))

  res=res[!is.na(res$iter),]

  if (!('name' %in% names(res))) 
    names(res)[1]='name'

  res=res[!is.na(res$name),]
  res=ddply(res, .(name,iter), function(x,log) data.frame(calcQ(x$stock,x$value)),log='log')
  its=max(as.numeric(ac(res$iter)))

  res.=transform(melt(res,id=c('name','iter')),params=paste(variable,name,sep=''))[,c('params','value','iter')]
  names(res.)[2]='data'

  #bug
  res=as(res.,'FLPar')[,1]
  #res=FLCore::iter(res,seq(its))
  units(res)='NA'
  res=res.[with(res.,order(iter,params)),]
  
  if (dims(object@params)$iter==1)
    object@params=propagate(object@params,its)
  t.=rbind(object@params,FLPar(as.FLQuant(res)[,1,drop=T]))
  dmns=dimnames(t.)
  names(dmns)=c('params','iter')
  t.=FLPar(array(t.,dim=unlist(lapply(dmns,length)),dimnames=dmns))
  units(t.)='NA'

  object@params=t.
  #object@params=FLPar(rbind(FLPar(object@params),FLPar(res)))
  
  object@params}

setMethod('setControl<-', signature(object='biodyn',value='FLQuant'), function(object,value,
  min=0.1,max=10.0) {
  
  setParams(object)<-value
  setControl(object,min=min,max=max)<-params(object)
  
  return(object)})


setMethod('setControl<-', signature(object='biodyn',value='FLQuants'), 
          function(object,value,min=0.1,max=10.0) {

  setParams(object)<-value
  setControl(object,min=min,max=max)<-params(object)
    
  return(object)})

setMethod('control<-', signature(object='biodyn',value='FLPar'), function(object,value) {
  object@control=value
    
  return(object)})

setMethod('control', signature(object='biodyn'), function(object) {
  object@control})

#' controlFn
#' @description A utility function to help set up the \code{control} slot in \code{biodyn} 
#' 
#' @param r a \code{numeric} value with best guess
#' @param k a \code{numeric} value with best guess      
#' @param p a \code{numeric} value with best guess, default=1      
#' @param b0 a \code{numeric} value with best guess default=1
#' @param phaseR a \code{numeric} value for phase, default=1
#' @param phaseK a \code{numeric}  value for phase, default=1,
#' @param phaseP a \code{numeric}  value for phase, default=-1,
#' @param phaseB0 a \code{numeric}  value for phase, default=-1,
#' @param min a \code{numeric} a multipler for the best guess  
#' @param max \code{numeric} a multipler for the best guess
#'
#' @return a \code{control} object
#'  
#' @export
#' @rdname controlFn
#' 
#' @seealso \code{\link{biodyn}}  \code{\link{control}} 
#'   
#' @examples
#' \dontrun{
#' simBiodyn()
#'    }
## utility function for setting control object
controlFn=function(r,       k,       p=1,      b0=1,
                   phaseR=1,phaseK=1,phaseP=-1,phaseB0=-1,
                   min=.5,  max=2){ 
  
  dmns=list(params=c('r','k','p','b0'),
            option=c('phase','min','val','max'),  
            iter  =1)
  
  res=FLPar(array(0,unlist(laply(dmns,length)),dmns))
  res[,'val'][]=c(r,k,p,b0)
  res[,'min']=res[,'val']*min
  res[,'max']=res[,'val']*max
  
  res[,'phase']=c(phaseR,phaseK,phaseP,phaseB0)
  
  res}  
  
#' priorFn
#' @description A utility function to help set up the \code{prior} slot in \code{biodyn}.
#'            
#' @param ... any other parameter
#' 
#' #' @return \code{list} with om, ...
#'  
#' @export

#' @rdname priorFn
#' 
#' @seealso \code{\link{biodyn}}
#' 
#' @examples
#' \dontrun{
#'    priorFn(r=c(weight=1,0.3,0.03)
#'    }
priorFn=function(...){
  
  args=list(...)
  if ('list' %in% is(args[[1]]))
    args=args[[1]]
  
  res=biodyn::biodyn()@priors
  
  for (i in dimnames(res)$params[dimnames(res)$params %in% names(args)]){
    if ('weight' %in% names(args[[i]])) res[i,c('weight')]=args[[i]]['weight']
    if ('a'      %in% names(args[[i]])) res[i,c('a')]     =args[[i]]['a']
    if ('b'      %in% names(args[[i]])) res[i,c('b')]     =args[[i]]['b']}
  
  res}  

