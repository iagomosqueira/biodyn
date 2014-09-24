#' randJack
#'
#' Simulates a \code{biodyn} object for a catch series, given the parameter estimates  in the \code{param} 
#' slot and variance covariance matrix 
#' http://young.physics.ucsc.edu/jackboot.pdf
#' 
#' @param n  \code{numeric} with number of simulations 
#' @param object \code{biodyn} 
#' @param ... other arguments
#' 
#' @return \code{biodyn} with estimates of stock based on catch time series
#' 
#' @export
#' @rdname randJack
#'
#' @aliases randJack-method randJack,numeric,biodyn-method
#' randJack, numeric,FLPar,FLPar-method
#' randJack, numeric,FLQuant,FLQuant-method
#' randJack,numeric,biodyn,biodyn-method
#'
#' @examples
#' \dontrun{
#' #simulate an object with known properties
#' bd=simBiodyn()
#' bd=window(bd,end=49)
#' 
#' #simulate a proxy for stock abundance
#' cpue=(stock(bd)[,-dims(bd)$year]+stock(bd)[,-1])/2
#' cpue=rlnorm(1,log(cpue),.2)
#' 
#' #set parameters
#' setParams(bd) =cpue
#' setControl(bd)=params(bd)
#' control(bd)[3:4,"phase"]=-1
#' 
#' #fit
#' bd=fit(bd,cpue)
#' 
#' bdJK=fit(bd,jackknife(cpue))
#' 
#' 
#' bdRnd  =randJack(100,bd,bdJK)
#' 
#' plot(randJack(100,stock(bd)[,40:45],stock(bdJK)[,40:45]))
#' }
setGeneric('randJack',   function(n,object,sim,...)    standardGeneric('randJack'))

randJackFn<-function(n,object,sim){
  res=jackSummary(params(object),params(sim))
  
  cov=res$cov[modelParams(model(object)),modelParams(model(object))]
  object@params=propagate(object@params,n)
  object@catch =propagate(object@catch, n)
  object@catch[]=object@catch[,,,,,1]
  object@stock =propagate(object@stock, n)
  
  nms=dimnames(cov)[[1]]
  nms=nms[aaply(cov,1,function(x) !all(is.na(x)))]
    
  #FLParBug in drop so need c()
  params(object)[nms]=
    t(maply(seq(n),function(x) 
         mvrnorm(1,c(params(object)[nms,x,1]),cov[nms,nms,drop=T])))
  
  rtn=fwd(object,catch=object@catch)
  
  return(rtn)}

setMethod('randJack', signature(n="numeric",object='biodyn',sim="biodyn"),  
           function(n,object,sim) randJackFn(n,object,sim))

setMethod('randJack', signature(n="numeric",object='FLQuant',sim="FLQuant"),  
          function(n,object,sim){ 
              res=jackSummary(object,sim)
               
              if (is.null(res$cov))
                rnorm(n,res$mean,res$se)
              else
                mvrnorm(n,res$hat,res$cov)
              })
          
setMethod('randJack', signature(n="numeric",object='FLPar',sim="FLPar"),  
          function(n,object,sim){ 
            res=jackSummary(object,sim)
            
            if (sum(dim(res$cov))>0)
              mvrnorm(n,res$hat,res$cov)
            else  
              rlnorm(n,res$hat,res$se)
          })
