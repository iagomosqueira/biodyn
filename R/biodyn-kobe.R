utils::globalVariables(c('ddply','.','year','pctl','cast','kobeP','sims'))

#' kobe
#' 
#' @description Creates time series of stock relative to BMSY and harvest rate relative
#' to FMSY
#' 
#' @name kobe
#' 
#' @param object biodyn object
#' @param method missing
#' @param ... other arguments
#' 
#' @return data.frame or list of data.frames
#' 
#' @aliases kobe kobe,biodyns,missing-method kobe,biodyns,ANY-method kobe,biodyn,missing-method 
#' 
#' @export
#' 
#' @rdname kobe
#' 
#' @examples 
#' \dontrun{
#' sim()
#' }

#if (!isGeneric('kobe')) 
setGeneric('kobe',  function(object,method,...) standardGeneric('kobe'))

setMethod('kobe', signature(object='biodyn',method='missing'),
          function(object,what=c('sims','trks','pts','smry','wrms')[1],probs=c(0.75,0.5,.25),
                   year=NULL,nwrms=10,sim=NULL){
            
            if (substr(sim,1,4)=="jack"){
              stock  =stock(  object)%/%bmsy(object)
              harvest=harvest(object)%/%fmsy(object)
            }else if ((dims(ple4)$iter)==1){
              stock  =stock(  object)%/%bmsy(object)
              harvest=harvest(object)%/%fmsy(object)
            }else{          
              res=model.frame(mcf(FLQuants(stock   =stock(  object)%/%bmsy(object),
                                           stock   =stock(  object)%/%bmsy(object))))
            }
            
            if ('pts' %in% what & is.null(year)) year=range(object)['maxyear']-1
            kobeFn(res,what,probs,year,nwrms)})

setMethod('kobe', signature(object='biodyns'),
          function(object,what=c('sims','trks','pts','smry','wrms')[1],probs=c(0.75,0.5,.25),year=NULL,nwrms=10){
            
            res=ldply(object,  function(x) model.frame(mcf(FLQuants(stock  =stock(  x)%/%bmsy(x),
                                                                    harvest=harvest(x)%/%fmsy(x)))))
            
            kobeFn(res,what,probs,year,nwrms)})

kobeMar=function(x,ds=seq(0,4,.001)){
  
  mar=rbind(cbind(what='stock',
                  data.frame(value  =ds,
                             density=dnorm(ds,x@mng['bbmsy','hat'],x@mng['bbmsy','sd']))),
            cbind(what='harvest',
                  data.frame(value  =ds,
                             density=dnorm(ds,x@mng['ffmsy','hat'],x@mng['ffmsy','sd']))))
  return(mar)}


kobeFn=function(object,what=c('sims','trks','pts','smry','wrms')[1],
                probs=c(0.75,0.5,.25),year=NULL,nwrms=10){         
  
  trks. =NULL
  pts.  =NULL
  smry. =NULL
  wrms. =NULL
  sims. =NULL
  
  ## trks
  if ('trks' %in% what){
    
    trks.=rbind(ddply(object,.(year), function(x) data.frame(quantity='stock',  pctl=probs,value=quantile(x$stock,    probs, na.rm=TRUE))),
                ddply(object,.(year), function(x) data.frame(quantity='harvest',pctl=probs,value=quantile(x$harvest,  probs, na.rm=TRUE))))
    
    trks.=transform(trks.,pctl=paste(substr(ac(signif(pctl,2)),3,nchar(ac(signif(pctl,2)))),ifelse(nchar(ac(trks.$pctl))==3,'0',''),'%',sep=''))
    trks.=cast(trks.,year+pctl~quantity,value='value') 
  }
  
  if ('pts' %in% what & !is.null(year))
    pts. =object[object$year==year,]
  
  if ('smry' %in% what)
    smry. =ddply(kobeP(sims), .(year), function(x) data.frame(stock      =median(stock(object),       na.rm=TRUE),
                                                              harvest    =median(harvest(object),     na.rm=TRUE),
                                                              red        =mean(  x$red,         na.rm=TRUE),
                                                              yellow     =mean(  x$yellow,      na.rm=TRUE),
                                                              green      =mean(  x$green,       na.rm=TRUE),
                                                              overFished =mean(  x$overFished,  na.rm=TRUE),
                                                              overFishing=mean(  x$overFishing, na.rm=TRUE)))
  if ('wrms' %in% what){          
    wrms =sample(unique(res$iter),nwrms)
    wrms.=sims[sims$iter %in% wrms,]
  }
  
  if ('sims' %in% what)     
    sims. =object
  
  res=list(trks=trks.,pts=pts.,smry=smry.,wrms=wrms.,sims=sims.)
  
  if (length(what)==1) res[[what]] else res[what]}