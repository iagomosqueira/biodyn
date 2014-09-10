#' rand
#'
#' Simulates a \code{biodyn} object for a catch series, given the parameter estimates  in the \code{param} 
#' slot and variance covariance matrix 
#' http://young.physics.ucsc.edu/jackboot.pdf
#' 
#' @param object of class \code{biodyn} 
#' @param catch of class \code{FLQuant} with the catch time series.
#' @param ... other arguments
#' 
#' @return \code{biodyn} with estimates of stock based on catch time series
#' 
#' @export
#' @rdname rand
#'
#' @aliases rand-method rand,biodyn,FLQuant-method
#'
#' @examples
#' \dontrun{
#' bd=simBiodyn()
#' cpue=(stock(bd)[,-dims(bd)$year]+stock(bd)[,-1])/2
#' cpue=rlnorm(1,log(cpue),.2)
#' bd  =fit(bd,jackknife(cpue))
#' bd  =jackRand(100,object)
#' }
setGeneric('jackRand',   function(n,object,...)    standardGeneric('jackRand'))

jackRand<-function(n,object){
  res=jackSummary(params(object))
  
  vcov=as.matrix(vcov(object)[,,1,drop=T])
  
  object     =iter(object,1)
  vcov[]=0
  diag(vcov)=res$se^2
  
  object@vcov=FLPar(vcov)
  catch=catch(object)
  cov  =vcov(object)[biodyn:::modelParams(model(object)),
                     biodyn:::modelParams(model(object))]

  object@params=propagate(object@params,n)
  object@catch =propagate(object@catch, n)
  object@catch[]=object@catch[,,,,,1]
  object@stock =propagate(object@stock, n)
  
  nms=dimnames(cov)[[1]]
  nms=nms[aaply(cov,1,function(x) !all(is.na(x)))]
    
  #FLParBug in drop so need c()
  params(object)[nms]=
    t(maply(seq(n),function(x) 
         mvrnorm(1,c(params(object)[nms,x,1]),cov[nms,nms,1,drop=T])))
  
  object=fwd(object,catch=catch)
  
  return(object)}

setMethod('rand', signature(n="numeric",object='biodyn'),  
          function(n,object) jackRand(n,object))
