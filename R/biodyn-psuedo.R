#' psuedo
#'
#' Simulates a \code{biodyn} object for a catch series, given the parameter estimates  in the \code{param} 
#' slot and variance covariance matrix 
#' 
#' @param object of class \code{biodyn} 
#' @param catch of class \code{FLQuant} with the catch time series.
#' @param ... other arguments
#' 
#' @return \code{biodyn} with estimates of stock based on catch time series
#' 
#' @export
#' @rdname psuedo
#'
#' @aliases psuedo-method psuedo,biodyn,FLQuant-method
#'
#' @examples
#' \dontrun{
#'    bd=psuedo(bd,FLQuant(10000,dimnames=list(year=2009:2012)))
#'    plot(bd)}
setGeneric('psuedo',   function(object,catch,...)    standardGeneric('psuedo'))

setMethod('psuedo', signature(object='biodyn',catch='FLQuant'),  
    function(object,catch,cov=vcov(object)[modelParams(model(object)),modelParams(model(object))],rand=TRUE){

  if (dims(catch)$iter>1 & dims(object)$iter==1)
     object=propagate(object,dims(catch)$iter)
  
  if (rand){
    nms=dimnames(cov)[[1]]
    nms=nms[aaply(cov,1,function(x) !all(is.na(x)))]
    
    #FLParBug in drop so need c()
    params(object)[nms]=t(maply(seq(dims(object)$iter),function(x) 
          mvrnorm(1,c(params(object)[nms,x]),cov[nms,nms,x,drop=T])))
    }
 
  object=fwd(object,catch=catch)

  return(object)})
