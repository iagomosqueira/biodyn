utils::globalVariables(c('ggplot','geom_line','aes','yield',
                         'geom_point','cast','xlab','ylab',
                         'jack.mean','facet_wrap','geom_ribbon',
                         'jack.se','expand_limits','jack.bias',
                         'theme','element_blank','y'))

utils::globalVariables(c('laply','ldply','geom_path','.id',
                         'scale_size_manual','scale_linetype_manual'))

utils::globalVariables(c('geom_path','scale_size_manual','scale_linetype_manual'))

utils::globalVariables(c('qname','What'))


whooow  =function(x,fn,probs)
  as.data.frame(FLQuants(lapply(fn,
                                function(fn,x)
                                  quantile(fn(x), probs=probs, type=type, na.rm=T), x=x))) 


#' plot
#' 
#' Creates a \code{ggplot2} object that plots time series of biomass, harvest rate and catch. The basic object can then be modified by adding ggplot2 layers.
#'
#' @aliases plot,biodyns,missing-method
#' @param x an object of class \code{biodyn} 
#' @param worm iters
#' @param y second argument
#' @param probs numeric vector of probabilities with values in [0,1]. (Values up to 2e-14 outside that range are accepted and moved to the nearby endpoint.)  
#' @param na.rm FALSE 
#' @param type an integer between 1 and 9 selecting one of the nine quantile algorithms detailed below to be used.
#' @param fn  functions
#' @param facet facet for panels
#'
#' @importFrom reshape cast
#' @import ggplot2
#' 
#' @return an \code{ggplot2} object
#' 
#' @seealso \code{\link{plotPrd}}, \code{\link{plotEql}}
#' 
#' @export
#' @rdname plot
#' 
#'
#' @examples
#' \dontrun{
#' 
#' x =simBiodyn()
#' x =window(x,end=49)
#' bd1=fwd(x,harvest=rlnorm(200,log(harvest(x)[,-1]),.2))
#' bd2=fwd(x,harvest=rlnorm(200,log(harvest(x)[,-1])*1.5,.2))
#' plot(biodyns("1"=bd1,"2"=bd2))
#' } 
setMethod('plot', signature(x='biodyn', y='missing'),
  function(x, y, probs=c(0.95,0.75,0.50,0.25,0.05), na.rm=FALSE, type = 7, 
    worm =NULL,
    fn   =list('Stock'  =function(x) stock(x), 
              'Harvest'=function(x) harvest(x),
              'Yield'  =function(x) catch(x)),
    facet=facet_wrap(~qname,scales='free',ncol=1),...){
    
      ## contains iters
      if (dims(x)$iter>=length(probs)){  
        res=whooow(x,fn,probs)
        res=transform(res,iter=factor(iter,labels= paste(as.integer(probs*100),"%",sep="")))
        res=cast(res,quant+year+unit+season+area+qname~iter,value="data")
        
        nprob=length(probs)
        if (any(0.5 %in% probs)){
          val=grep("%",names(res))

          mdn=res[,c(1:6,val[ (val==grep("50%",names(res)))])]
          
          res=res[,c(1:6,val[!(val==grep("50%",names(res)))])]}
        
        rbn=mdply(data.frame(ymax=dim(res)[2]-seq(nprob/2)+1,
                               ymin=6+seq(nprob/2)),
                    function(ymax,ymin){
                       data.frame(res[,c(1:6)],ymin=res[,ymin],ymax=res[,ymax],
                                  CI=paste(names(res)[ymin],names(res)[ymax],sep="-"))
                    })
          
          p=ggplot(rbn)+
            geom_ribbon(aes(year,ymax=ymax,ymin=ymin,alpha=CI),fill="red")
        
         if (any(0.5 %in% probs))
           p=p+geom_line(  aes(year,`50%`),data=mdn)

      }else{
        res=whooow(x,fn,0.5)
        p  =ggplot(res) + geom_path(aes(x=year,y=data))    }
      
      p = p + facet + 
            expand_limits(y = 0) +
            xlab('Year') + ylab('')
      
      if (length(worm) > 0)
        if (length(worm)<=dims(x)$iter)
          if (!(length(worm)==1 & is.na(worm[1])))  
            p=p+geom_path(aes(year,data,group=iter,colour=iter),
                            data=transform(subset(as.data.frame(FLQuants(lapply(fn,function(f,x) f(x), x=x))),iter %in% worm),iter=factor(iter)))
      
      p}) 
    
setMethod('plot', signature(x='biodyns', y='missing'),
  function(x, y, probs=c(0.95,0.75,0.50,0.25,0.05),type=7,na.rm=FALSE,
    facet=facet_wrap(~qname,scales='free',ncol=1),
           fn=list('Stock'  =function(x) stock(x), 
                   'Harvest'=function(x) harvest(x),
                   'Yield'  =function(x) catch(x)),...)
    {  
    res=ldply(x,function(x) plot(x,probs=probs,type=type,na.rm=na.rm)$data)
    mdn=ldply(x,function(x) whooow(x,fn,probs=.5))
    
    if (length(unique(res$CI))>=1){
      p=ggplot(res)+
        geom_ribbon(aes(year,ymax=ymax,ymin=ymin,alpha=CI,group=paste(.id,CI),fill=.id))+
        geom_line(  aes(year,data,col=.id,group=.id),data=mdn)+
        scale_alpha_manual(values=seq(.1,.75,length.out=length(unique(res$CI))))
      
    }else{
      p=ggplot(res)+
        geom_line(aes(year,data,col=.id,group=.id),data=mdn)
      
    }
    
    p = p + facet + 
      expand_limits(y = 0) +
      xlab('Year') + ylab('')
  
    p})

# @param  \code{fn}, a list of functions that estimate the quantities for plotting
# @param  \code{probs}, a vector specifying the percentiles for plotting, these are c(0.95,0.50,0.05) by default.
# @param  \code{size}, thinkness of percentile lines
# @param  \code{lty}, line type for percentiles
# @param \code{facet}, a layer that determines the facetting of the plot


#' plotPrd
#' 
#' Creates a \code{ggplot2} object that plots equilibrium values of biomass, harvest rate and catch against each other.
#' The basic object can then be modified by adding ggpot2 layers.
#'
#' @param  x an object of class \code{biodyn} 
#' @param  biomass optional argument, an FLQuant with biomass at beginning of year 
#' @param ... other arguments
#'
#' @return an \code{ggplot2} object
#' 
#' @seealso \code{\link{plotMSE}}, \code{\link{plotEql}}
#' 
#' @export
#' @rdname plotPrd
#'
#' @aliases plotPrd,biodyn,FLBRP-method plotPrd,biodyn,FLQuant-method plotPrd,biodyn,missing-method
#' 
#' @examples
#' \dontrun{
#'  refpts('logistic',FLPar(msy=100,k=500))
#' }
setGeneric('plotPrd', function(x,biomass,...) standardGeneric('plotPrd'))
setMethod('plotPrd',signature(x='biodyn',biomass='missing'),  
          function(x,biomass,...) plotPrdfn(x,biomass=FLQuant(seq(0,max(params(x)['k']),length.out=101)),...))
setMethod('plotPrd',signature(x='biodyn',biomass='FLQuant'),  
          function(x,biomass,...) plotPrdfn(x,biomass,...))
# setMethod('plotPrd',signature(x='biodyn',biomass='FLBRP'),  
#           function(x,biomass,II=FALSE,...) plotProdfn(bd=x,brp=biomass,II=II,...))

plotPrdfn=function(x,biomass=FLQuant(seq(0,max(params(x)['k']),length.out=101)),...) {
  if ((dims(x)$iter>1 | dims(params(x))$iter>1) & dims(biomass)$iter==1) 
    biomass=propagate(biomass,max(dims(x)$iter,dims(params(x))$iter))
  
  p <-  ggplot(model.frame(FLQuants(stock=biomass, yield=FLQuant(computePrd(x,biomass))))) +
    geom_line(aes(stock, yield, group=iter, col=iter)) +
    geom_point(aes(bmsy,msy,col=iter),size=2,data=cast(as.data.frame(refpts(x)),iter~refpts,value='data')) +
    xlab('Stock') + ylab('Surplus Production')
  
  p} 

plotProdfn=function(bd,brp,II=FALSE){
  res=cbind(What='Age I',  model.frame(FLQuants('catch'=catch(brp),'stock'=stock(brp)),drop=T))
  
  if (II){
    landings.wt(brp)=stock.wt(brp)*mat(brp)
    res=rbind(res, cbind(What='Age II',   model.frame(FLQuants('catch'=catch(brp),'stock'=stock(brp)),drop=T)))}
  
  if (!is.null(bd)){
    bm =FLQuant(seq(0, max(params(bd)['k']), length.out = 101))
    res=rbind(res,cbind(What='Biomass', 
                        model.frame(mcf(FLQuants(catch=computePrd(bd,bm),stock=bm)),drop=T)))}
  
  ggplot(res)}

#' plotEql
#' 
#' Creates a \code{ggplot2} object that plots time series of biomass, harvest rate and catch. The basic object can then be modified by adding ggpot2 layers.
#'
#' @param  x an object of class \code{biodyn} 
#' @param  biomass an object of holding biomass at beginning of year 
#' @param ... other arguments
#'
#' @return an \code{ggplot2} object
#' 
#' @seealso \code{\link{plotPrd}}\code{\link{plotMSE}}
#' 
#' @export
#' @rdname plotEql
#'
#' @aliases plotEql,biodyn,FLBRP-method plotEql,biodyn,FLQuant-method plotEql,biodyn,missing-method
#' 
#' @examples
#' \dontrun{
#'  refpts('logistic',FLPar(msy=100,k=500))
#' }
setGeneric('plotEql', function(x,biomass,...) standardGeneric('plotEql'))
setMethod('plotEql',signature(x='biodyn',biomass='missing'),  
          function(x,biomass,...) plotEqlfn(x,biomass=FLQuant(seq(0,max(params(x)['k']),length.out=101)),...))
setMethod('plotEql',signature(x='biodyn',biomass='FLQuant'),  
          function(x,biomass,...) plotEqlfn(x,biomass,...))
# setMethod('plotEql',signature(x='biodyn',biomass='FLBRP'),  
#           function(x,biomass,II=FALSE,...) plotProdfn(bd=x,brp=biomass,II=II,...))

plotEqlfn=function(x,biomass=FLQuant(seq(0,max(params(x)['k']),length.out=101)),...) {
  if ((dims(x)$iter>1 | dims(params(x))$iter>1) & dims(biomass)$iter==1) 
    biomass=propagate(biomass,max(dims(x)$iter,dims(params(x))$iter))
  
  res=model.frame(FLQuants(stock=biomass, 
                           yield=FLQuant(computePrd(x,biomass))))
  res=transform(res,harvest=yield/stock)
  
  res=rbind(cbind(pnl="Equilibrium Stock v Harvest Rate",
                  transform(res,x=harvest,y=stock)[,-(7:9)]),
            cbind(pnl="Equilibrium Production v Harvest Rate",
                  transform(res,x=harvest,y=yield)[,-(7:9)]),
            cbind(pnl="Equilibrium Production v Stock",
                  transform(res,x=stock,  y=yield)[,-(7:9)]))
  
    p=ggplot(res)+geom_line(aes(x,y))+
      facet_wrap(~pnl,scales="free",ncol=1)+    
      xlab('') + ylab('')
  
  p} 

plotProdfn=function(bd,brp,II=FALSE){
  res=cbind(What='Age I',  model.frame(FLQuants('catch'=catch(brp),'stock'=stock(brp)),drop=T))
  
  if (II){
    landings.wt(brp)=stock.wt(brp)*mat(brp)
    res=rbind(res, cbind(What='Age II',   model.frame(FLQuants('catch'=catch(brp),'stock'=stock(brp)),drop=T)))}
  
  if (!is.null(bd)){
    bm =FLQuant(seq(0, max(params(bd)['k']), length.out = 101))
    res=rbind(res,cbind(What='Biomass', 
                        model.frame(mcf(FLQuants(catch=computePrd(bd,bm),stock=bm)),drop=T)))}
  
  ggplot(res)}


#' plotMSE
#'
#' Creates a \code{ggplot2} object that plots absolute and relative to MSY benchmarks time series of 
#' ssb, biomass, harvest rate and catch for  \code{FLStock} and  \code{biodyn}  objects
#' The basic object can then be modified by adding ggpot2 layers.
#'
#' @param  x, an object of class \code{biodyn} 
#' @param  y, an object of class \code{FLStock} 
#' @param  z, an object of class \code{FLBRP} 
#' @param ... other arguments
#'
#' @return an \code{ggplot2} object
#' 
#' @seealso \code{\link{plotPrd}}
#' 
#' @export
#' @rdname plotMSE
#'
#' @aliases plotMSE,biodyn,FLStock,FLBRP-method
#' 
#' @examples
#' \dontrun{ 
#' refpts('logistic',FLPar(msy=100,k=500))
#' }
setGeneric('plotMSE', function(x,y,z,...) standardGeneric('plotMSE'))
# setMethod('plotMSE',signature(x='biodyn',y='FLStock',z='FLBRP'),  
#           function(x,y,z,...) plotMSEfn(x,y,z,...))


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
    facet_wrap(qname~Type,scales='free',ncol=2)  }

##############################################################
#' plotJack
#' 
#' Create a \code{ggplot2} plot based on a jack knifed biodyn and plots 
#' time series of biomass and harvest rate. 
#' The basic object can then be modified by adding ggpot2 layers.
#'
#' @param  x an object of class \code{biodyn} that has been jack knifed, i.e. by 
#' providing a jack knifed CPUE series to fit
#' @param y the original \code{biodyn} object
#' @param ncol number of colums in plot panel
#'
#' @return an \code{ggplot2} object
#' 
#' @export
#' @rdname plotJack
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
#' plotJack(bdJK)
#' bd  =randJack(100,bd)
#' }
plotJack=function(x,y,ncol=1){
   
  df=rbind(cbind(qname='stock',  
                 model.frame(jackSummary(stock(  x),stock(  y)),drop=TRUE)),
           cbind(qname='harvest',
                 model.frame(mcf(jackSummary(harvest(x),harvest(y))),drop=TRUE)))
  
  # basic plot data vs. year
  p=ggplot(data=df, aes(x=year, y=mean))+
    facet_wrap(~qname,ncol=ncol,scales='free_y')+
    geom_ribbon(aes(x=year, 
                  ymin=mean-2*se, 
                  ymax=mean+2*se),
                fill='blue', alpha = .20)+
    # line + xlab + ylab + limits to include 0 +
    geom_line(colour='red') + xlab('Year') + ylab('') + expand_limits(y=0) +
    #geom_line(aes(year,jack.mean+jack.bias),colour='black')   +
    # no legend
    theme(legend.title = element_blank())
  
  return(p)}

