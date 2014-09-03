utils::globalVariables(c('ggplot','geom_line','aes','yield',
                         'geom_point','cast','xlab','ylab',
                         'jack.mean','facet_wrap','geom_ribbon',
                         'jack.se','expand_limits','jack.bias',
                         'theme','element_blank'))

utils::globalVariables(c('laply','ldply','geom_path','.id',
                         'scale_size_manual','scale_linetype_manual'))

utils::globalVariables(c('geom_path','scale_size_manual','scale_linetype_manual'))

utils::globalVariables(c('qname','What'))


##############################################################
#' plot
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
#' @aliases plot,biodyns,missing
#'
#' @examples
#' \dontrun{
#' refpts('logistic',FLPar(msy=100,k=500))
#' } 
setMethod('plot', signature(x='biodyn', y='missing'),
  function(x, y, probs=c(0.95,0.50,0.05), size=c(0.5,1.0,0.5), lty=c(2,1,2), worm=NA,
    facet=facet_wrap(~qname,scale='free',ncol=1),
    fn=list('Stock'  =function(x) stock(x), 
            'Harvest'=function(x) harvest(x),
            'Yield'  =function(x) catch(x)),...){

    plotComp(x,fn,probs,size,lty,facet,worm=worm)})

setMethod('plot', signature(x='biodyns', y='missing'),
  function(x, y, probs=c(0.95,0.50,0.05), size=c(0.5,1.0,0.5), lty=c(2,1,2),
    facet=facet_wrap(~qname,scale='free',ncol=1),
           fn=list('Stock'  =function(x) stock(x), 
                   'Harvest'=function(x) harvest(x),
                   'Yield'  =function(x) catch(x)),...)
    
    plotComps(x,fn,probs,size,lty,facet))

# @param  \code{fn}, a list of functions that estimate the quantities for plotting
# @param  \code{probs}, a vector specifying the percentiles for plotting, these are c(0.95,0.50,0.05) by default.
# @param  \code{size}, thinkness of percentile lines
# @param  \code{lty}, line type for percentiles
# @param \code{facet}, a layer that determines the facetting of the plot


##############################################################
#' Create a \code{ggplot2} plot
#'
#' Creates a \code{ggplot2} object that plots absolute and relative to MSY benchmarks time series of 
#' ssb, biomass, harvest rate and catch for  \code{FLStock} and  \code{biodyn}  objects
#' The basic object can then be modified by adding ggpot2 layers.
#'
#' @param  x, an object of class \code{biodyn} 
#' @param  y, an object of class \code{FLStock} 
#' @param  z, an object of class \code{FLBRP} 
#'
#' @return an \code{ggplot2} object
#' 
#' @seealso \code{\link{plotSP}}
#' 
#' @export
#' @docType methods
#' @rdname plotMSE
#'
#' @aliases plotMSE,biodyn,FLStock,FLBRP-method
#' 
#' @examples
#' \dontrun{ 
#' refpts('logistic',FLPar(msy=100,k=500))
#' }
setGeneric('plotMSE', function(x,y,z,...) standardGeneric('plotMSE'))
setMethod('plotMSE',signature(x='biodyn',y='FLStock',z='FLBRP'),  
          function(x,y,z,...) plotMSEfn(x,y,z,...))


## compares age and biomass based time series      
plotMSEfn=function(mp,om,brp){
  ### OM ######################################
  ## absolute
  omAbs=cbind(Type='Absolute',
              rbind(as.data.frame(FLQuants(om,'stock','ssb','catch'),          drop=T),
                    as.data.frame(FLQuants(harvest=catch(om)/stock(om)),drop=T)))
  
  ## relative
  omRel=cbind(Type='Relative',as.data.frame(mcf(
    FLQuants('stock'  =stock(om)%/%refpts(brp)['msy','biomass'],
             'ssb'    =ssb(  om)%/%refpts(brp)['msy','biomass'],
             'harvest'=fbar( om)%/%refpts(brp)['msy','harvest'],
             #'harvest'=catch(om)/stock(om)%/%(refpts(brp)['msy','yield']/refpts(brp)['msy','biomass']),
             'catch'  =catch(om)%/%refpts(brp)['msy','yield'])),drop=T))
  
  ### MP ######################################
  ## absolute
  mpAbs=cbind(Type='Absolute',as.data.frame(FLQuants(mp,'stock','harvest','catch'),drop=T))
  
  ## relative
  mpRel=cbind(Type='Relative',as.data.frame(FLQuants('stock'  =stock(  mp)%/%bmsy(mp),
                                                     'harvest'=harvest(mp)%/%fmsy(mp),
                                                     'catch'  =catch(  mp)%/%msy( mp)),drop=T))
  
  ggplot(subset(rbind(cbind(What='OM',rbind(omAbs,omRel)),
                      cbind(What='MP',rbind(mpAbs,mpRel))),qname!='ssb'))+
    geom_line(aes(year,data,group=What,col=What))+
    facet_wrap(qname~Type,scale='free',ncol=2)  }


##############################################################
#' plotSP
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
#' @aliases plotSP,biodyn,FLBRP-method plotSP,biodyn,FLQuant-method plotSP,biodyn,missing-method
#' 
#' @examples
#' \dontrun{
#'  refpts('logistic',FLPar(msy=100,k=500))
#' }
setGeneric('plotSP', function(object,biomass,...) standardGeneric('plotSP'))
setMethod('plotSP',signature(object='biodyn',biomass='missing'),  
          function(object,biomass,...) plotSPfn(object,biomass=FLQuant(seq(0,max(params(object)['k']),length.out=101)),...))
setMethod('plotSP',signature(object='biodyn',biomass='FLQuant'),  
          function(object,biomass,...) plotSPfn(object,biomass,...))
setMethod('plotSP',signature(object='biodyn',biomass='FLBRP'),  
          function(object,biomass,II=FALSE,...) plotProdfn(bd=object,brp=biomass,II=II,...))

plotSPfn=function(object,biomass=FLQuant(seq(0,max(params(object)['k']),length.out=101)),...) {
  if ((dims(object)$iter>1 | dims(params(object))$iter>1) & dims(biomass)$iter==1) 
    biomass=propagate(biomass,max(dims(object)$iter,dims(params(object))$iter))
  
  p <-  ggplot(model.frame(FLQuants(stock=biomass, yield=FLQuant(computeSP(object,biomass))))) +
    geom_line(aes(stock, yield, group=iter, col=iter)) +
    geom_point(aes(bmsy,msy,col=iter),size=2,data=cast(as.data.frame(refpts(object)),iter~refpts,value='data')) +
    xlab('Stock') + ylab('Surplus Production')
  print(p)
  invisible(p)} 

plotProdfn=function(bd,brp,II=FALSE){
  res=cbind(What='Age I',  model.frame(FLQuants('catch'=catch(brp),'stock'=stock(brp)),drop=T))
  
  if (II){
    landings.wt(brp)=stock.wt(brp)*mat(brp)
    res=rbind(res, cbind(What='Age II',   model.frame(FLQuants('catch'=catch(brp),'stock'=stock(brp)),drop=T)))}
  
  if (!is.null(bd)){
    bm =FLQuant(seq(0, max(params(bd)['k']), length.out = 101))
    res=rbind(res,cbind(What='Biomass', 
                        model.frame(mcf(FLQuants(catch=computeSP(bd,bm),stock=bm)),drop=T)))}
  
  ggplot(res)}

##############################################################
#' plotJack
#' 
#' Create a \code{ggplot2} plot based on a jack knifed biodyn and plots 
#' time series of biomass and harvest rate. 
#' The basic object can then be modified by adding ggpot2 layers.
#'
#' @param  \code{x}, an object of class \code{biodyn} that has been jack knifed, i.e. by 
#' providing a jack knifed CPUE series to fit
#' @param \code{y} the original \code{biodyn} object
#'
#' @return an \code{ggplot2} object
#' 
#' @export
#' @docType methods
#' @rdname plotJack
#'
#' @examples
#' \dontrun{
#' bd=fit(bd,cpue)
#' jk=fit(bd,jacknife(cpue))
#' plotJack(jk,bd)
#' }  
plotJack=function(x,y,ncol=1){
  
  js=function(x,y, ...) {
    
    n <- dims(x)$iter 
    
    mnU <- apply(x, 1:5, mean)   
    
    SS <- apply(sweep(x, 1:5, mnU,'-')^2, 1:5, sum)
    
    bias <- (n - 1) * (mnU - y)
    se <- sqrt(((n-1)/n)*SS)
    
    res=FLQuants(list(jack.mean=y, jack.se=se, jack.bias=bias))
    
    attributes(res)$jackknife=TRUE
    
    return(res)}
  
  df=rbind(cbind(qname='stock',  
                 model.frame(js(stock(  x),stock(  y)),drop=TRUE)),
           cbind(qname='harvest',
                 model.frame(js(harvest(x),harvest(y)),drop=TRUE)))
  
  # basic plot data vs. year
  p=ggplot(data=df, aes(x=year, y=jack.mean))+
    facet_wrap(~qname,ncol=ncol,scale='free_y')+
    geom_ribbon(aes(x=year, ymin=jack.mean-2*jack.se, 
                    ymax=jack.mean+2*jack.se),
                fill='blue', alpha = .20)+
    # line + xlab + ylab + limits to include 0 +
    geom_line(colour='red') + xlab('Year') + ylab('') + expand_limits(y=0) +
    geom_line(aes(year,jack.mean+jack.bias),colour='black')   +
    # no legend
    theme(legend.title = element_blank())
  
  return(p)}

