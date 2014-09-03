##############################################################
#' power
#'
#' Estimates parameters in a \code{biodyn} class by fitting catch to CPUE indeices
#' 
#' @param   \code{object}, an object of class \code{biodyn}
#' @param   \code{test}, a \code{data.frame} with 3 columns \code{param,q,lower.tail} for the choice of parameter, the desired qunatile and teh tail of the distribution for the test.
#'
#' @return an \code{FLPar} object with probability of being greater (lower.tail=FALSE) or less (lower.tail=TRUE) then q
#'
#' @aliases powerAnalysis-method powerAnalysis,biodyn,biodyn-method  powerAnalysis,biodyn,missing-method
#' 
#' @export
#' @docType methods
#' @rdname power
#'
setGeneric('powerAnalysis',  function(object,ref,...) signature(object,ref)) 
          
setMethod('powerAnalysis',  signature(object='biodyn',ref='missing'), 
          function(object,ref=ref,test=data.frame(param =c('bbmsy','bk','bratio','ffmsy', 'fr','fratio','slopeb','slopef'),
                                              q         =c(     1,    0,       1,      0,    1,       0,       1,      0),
                                              lower.tail=c( FALSE, TRUE,   FALSE,   TRUE,FALSE,    TRUE,   FALSE,    TRUE))){
  
  its=dims(object)$iter
  
  if (its==1) object=propagate(object,2)
  
  res=maply(seq(dim(test)[1]),function(x,mng,test) {
    par=ac(test[x,'param'])
    pnorm(test[x,'q'],mng[par,'hat'],mng[par,'sd'],test[x,'lower.tail'])
    
  },mng=object@mng,test=test)
  res=as(res,'FLPar')
  dimnames(res)$X1=test$param
  names(dimnames(res))='param'
  res=as(res,'FLPar')
  
  if (its==1) res=FLCore::iter(object,1)
  
  return(res)})

setMethod('powerAnalysis',  signature(object='biodyn',ref='biodyn'), 
          function(object,ref,test=data.frame(param     =c('bnow','fnow','bnow','fnow'),
                                              ref       =c('bmsy','fmsy','bnow','fnow'),
                                              q         =c(     1,    0,      1,    0),
                                              lower.tail=c( FALSE, TRUE, FALSE, TRUE))){
  
  its=dims(object)$iter
            
  if (its==1) {
    object=propagate(object,2)
    ref   =propagate(ref,   2)}
    
  if (dims(object)$iter==1 & dims(ref)$iter>1)
    object=propagate(object,dims(ref)$iter)
  if (dims(object)$iter>1 & dims(ref)$iter==1)
    ref=propagate(ref,dims(object)$iter)
  
  
  obj=object@mng
  ref=ref@mng
  sd =obj[,'sd']^2*(1/ref[,'hat'])^2 + ref[,'sd']^2*(obj[,'hat']/ref[,'hat']^2)^2
  hat=obj[,'hat']/ref[,'hat']
  mng=object@mng
  mng[,'sd'] =sd
  mng[,'hat']=hat
  
  res=maply(seq(dim(test)[1]),function(x,mng,test) {
    par=ac(test[x,'param'])
    pnorm(test[x,'q'],mng[par,'hat'],mng[par,'sd'],test[x,'lower.tail'])
    
  },mng=mng,test=test)
  
  as(res,'FLPar')
  dimnames(res)$X1=test$param
  names(dimnames(res))='param'
  res=as(res,'FLPar')
  
  if (its==1){
    object=FLCore::iter(object,1)
    ref   =FLCore::iter(ref,   1)}
    
  return(res)})
