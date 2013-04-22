utils::globalVariables(c("ggplot","geom_line","aes","yield","geom_point","cast","xlab","ylab"))

##############################################################
#' Create a \code{ggplot} plot
#'
#' Creates a \code{ggplot2} object that plots time series of biomass, harvest rate and catch. The basic object can then be modified by adding ggpot2 layers.
#'
#' @param  \code{x}, an object of class \code{biodyn} 
#'
#' @return an \code{ggplot2} object
#' 
#' @seealso \code{\link{plotSP}} 
#' 
#' @export
#' @docType methods
#' @rdname plot
#'
#' @examples
#' refpts("logistic",FLPar(msy=100,k=500))
#'  
setMethod("plot", signature(x="biodyn", y="missing"),
  function(x, y, probs=c(0.95,0.50,0.05), size=c(0.5,1.0,0.5), lty=c(2,1,2), worm=NA,
    facet=facet_wrap(~qname,scale="free",ncol=1),
    fn=list("Stock"=stock, "Harvest"=function(x) catch(x)/stock(x)[,dimnames(catch(x))$year],"Yield"=function(x) catch(x)),...){

    plotComp(x,fn,probs,size,lty,facet,worm=worm)})

setMethod("plot", signature(x="biodyns", y="missing"),
  function(x, y, probs=c(0.95,0.50,0.05), size=c(0.5,1.0,0.5), lty=c(2,1,2),
    facet=facet_wrap(~qname,scale="free",ncol=1),
    fn=list("Stock"=stock, "Harvest"=function(x) catch(x)/stock(x)[,dimnames(catch(x))$year],"Yield"=catch),...)
   
    plotComps(x,fn,probs,size,lty,facet))

# @param  \code{fn}, a list of functions that estimate the quantities for plotting
# @param  \code{probs}, a vector specifying the percentiles for plotting, these are c(0.95,0.50,0.05) by default.
# @param  \code{size}, thinkness of percentile lines
# @param  \code{lty}, line type for percentiles
# @param \code{facet}, a layer that determines the facetting of the plot

##############################################################
#' Create a \code{ggplot} plot
#'
#' Creates a \code{ggplot2} object that plots time series of biomass, harvest rate and catch. The basic object can then be modified by adding ggpot2 layers.
#'
#' @param  object, an object of class \code{biodyn} 
#' @param  biomass, an object of holding biomass at beginning of year 
#'
#' @return an \code{ggplot2} object
#' 
#' @seealso \code{\link{plotSP}}
#' 
#' @export
#' @docType methods
#' @rdname plotSP
#'
#' @examples
#' refpts("logistic",FLPar(msy=100,k=500))
#' 
plotSP=function(object,biomass=FLQuant(seq(0,max(params(object)["k"]),length.out=101))) {
   if ((dims(object)$iter>1 | dims(params(object))$iter>1) & dims(biomass)$iter==1) 
      biomass=propagate(biomass,max(dims(object)$iter,dims(params(object))$iter))

   p <-  ggplot(model.frame(FLQuants(stock=biomass, yield=FLQuant(computeSP(object,biomass))))) +
             geom_line(aes(stock, yield, group=iter, col=iter)) +
             geom_point(aes(bmsy,msy,col=iter),size=2,data=cast(as.data.frame(refpts(object)),iter~refpts,value="data")) +
             xlab("Stock") + ylab("Surplus Production")
   print(p)
   invisible(p)} 
