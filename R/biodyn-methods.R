setMethod('harvest', signature(object='biodyn'),
             function(object,when=.5,...) {
             
             yrs1=  dimnames(stock(object))$year
             yrs2=c(dimnames(stock(object))$year[-1],as.numeric(max(dimnames(stock(object))$year))+1)
             
             #res <- catch(object)/(stock(object)[,yrs1]*(1-when)+
             #                        stock(object)[,yrs2]*when)
             
             yrs=dimnames(catch(object))$year[dimnames(catch(object))$year %in% dimnames(catch(object))$year]
             res <- catch(object)[,yrs]/stock(object)[,yrs]
             units(res) <- "hr"
             return(res)
           })

setMethod('stock', signature(object='biodyn'),
          function(object,when=0) {
            
            when=max(min(when,1),0)
            if (when<=0) return(object@stock)
            
            yrs =  dimnames(stock(object))$year
            yrs1=  rev(rev(yrs)[-1])
            yrs2=  yrs[-1]
             
            (1-when)*stock(object)[,yrs1]+when*stock(object)[,yrs2]})

##############################################################
#' Calculates surplus production
#'
#' Calculates the surplus production for a biomass dynamic model given a level of stock biomass
#' 
#' @param  \code{object}, an object of class \code{biodyn} 
#'
#' @param \code{biomass}, stock biomaas, may be a \code{numerix},  \code{FLQuant} or missing. In the latte case the stock slot will be used.
#'
#' @return an \code{FLPar} object
#' 
#' @seealso \code{\link{plotSP}}
#' 
#' @export
#' @docType methods
#' @rdname sp
#'
#' @examples
#' \dontrun{ computeSP(bd,seq(0,params(bd)["k"])) }
#'  
setGeneric('computeSP',function(object,biomass,...) standardGeneric('computeSP'))
setMethod( 'computeSP', signature(object="biodyn",   biomass="missing"),     function(object,biomass=stock(object))  biodyn:::spFn(model(object),params(object),biomass))
setMethod( 'computeSP', signature(object="biodyn",   biomass="numeric"),     function(object,biomass)                biodyn:::spFn(model(object),params(object),biomass))
setMethod( 'computeSP', signature(object="biodyn",   biomass="FLQuant"),     function(object,biomass)                biodyn:::spFn(model(object),params(object),biomass))


# calcLogLik

calcSigma <- function(obs,hat=rep(0,length(obs)),error="log"){
  yrs=dimnames(obs)$year
  yrs=yrs[yrs %in% dimnames(hat)$year]
  hat=hat[,yrs]
  obs=obs[,yrs]
  
  if (error=="log"){
    hat=log(hat)
    obs=log(obs)}
  
  SS =sum((obs-hat)^2,na.rm=T)
  
  return((SS/length(hat))^.5)}

loglFn<-function(obs,se,hat=rep(0,length(obs))){
  flag=!is.na(obs) & !is.na(hat)
  obs =obs[flag]
  hat =hat[flag]
  
  SS<-sum((obs-hat)^2)
  
  n   <-length(obs)
  res <-(log(1/(2*pi))-n*log(se)-SS/(2*se^2))/2
  
  return(res)}

calcLogLik<-function(obs,hat=rep(0,length(obs)),error="log",type=1){
  
  yrs=dimnames(obs)$year
  yrs=yrs[yrs %in% dimnames(hat)$year]
  hat=hat[,yrs]
  obs=obs[,yrs]
  
  if (error=="log"){
    hat=log(hat)
    obs=log(obs)}
  
  se<-calcSigma(obs,hat)
  
  if (type==1) return(loglFn(se,obs,hat)) else
    if (type==2) return(-sum(dnorm(obs, hat, se, log=(error=="log"), na.rm=TRUE))) else
      if (type==3) return(sum((obs-hat)^2))}
