#' mvn
#' 
#' @description 
#' Takes a fitted biodyn object and uses the covariance matrix (in the \code{vcov} slot)  and the
#' parameter estimates (\code{params} slot) to use Monte Carlo simulation to generate new 
#' parameters
#'      
#' @param object, a \code{biodyn} object
#' @param n, \code{numeric} with number of \code{iter} to create
#' @param fwd, \code{logical} do you want to simuate historic time series of stock biomass?, default is \code{FALSE} 
#' @return \code{biodyn} with simuated time series 
#' @export
#' @docType functions
#' @rdname mvn
#' 
#' @examples
#' x=1
mvn=function(object,n,nms=dimnames(object@control[object@control[,"phase",]>0,])$params,
                      fwd=FALSE){
  
  res=mvrnorm(n,params(object)[nms,drop=T],vcov(object)[nms,nms,drop=T])
  dmns=list(iter=seq(n),params=nms)
  res=array(res,dim=unlist(llply(dmns,length)),dimnames=dmns)
  res=FLPar(t(res))
  units(res)="NA"
  
  params(object)=propagate(params(object),n)
  stock(object) =propagate(stock( object),n)
  
  params(object)[nms]=res
  
  if (fwd) object=fwd(object,catch=catch(object))
              
  object}
