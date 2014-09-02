utils::globalVariables(c('profileGrid'))
utils::globalVariables('maply')

#' profile
#'
#' Profiles biodyn  
#'
#' @param   \code{fitted}, an object of class \code{biodyn}
#' @param   \code{index}, an \code{FLQuant}, \code{FLQuants} or  \code{data.frame} object with CPUE indices
#' @param   \code{cmdOps}, a character string giving ADMB options see \url{http://www.admb-project.org/documentation/manuals/ADMBrefcard-A4.pdf/view}
#'
#' @export
#' @docType methods
#' @rdname fit
#'
#' @examples
#' \dontrun{
#' library(aspic)
#' library(biodyn)
#' data(asp)
#' }

# ### debugging stuff
# data(bd)
# fitted=biodyn(factor('pellat'),params(bd),catch=catch(bd))
# index=rlnorm(1,log(stock(bd)),.2)[,-60]
# setParams(fitted)     =index
# 
# 
# attach(list(maxsteps=11, range=0.5, ci=c(0.25, 0.5, 0.75, 0.95),
#             plot=TRUE,fixed=c()))
# which='r'
# fixed=c('p','b0')
# ###
# res=profile(bd,which='r',fixed=c('b0','p'),index,range=c(1.2,3.0))
# ggplot(res)+geom_line(aes(r,ll))
# v <- ggplot(res, aes(r, k, z = ll))
# v <- ggplot(res, aes(r, k, z = ll))+ stat_contour(aes(colour = ..level..), size = 1)
setMethod('profile', signature(fitted='biodyn'),
      function(fitted,index,
                   which,
                   range=seq(0.5,1.5,length.out=21),
                   fn   =function(x) cbind(model.frame(params(x)),ll=model.frame(x@ll)[,1],model.frame(refpts(x))[,-4],
                                           stock  =c(stock(  x)[,ac(range(x)['maxyear'])]%/%bmsy(x)),
                                           harvest=c(harvest(x)[,ac(range(x)['maxyear'])]%/%fmsy(x))),
                   run  =TRUE,...){
  
        if (dims(index)$maxyear>=dims(stock(fitted))$maxyear) stop('index years greater in length than stock')
        if (dims(fitted)$iter>1) stop('can only be done for a single iter')
                 
        if (dim(fitted@control)[3]==1){
           fitted@control=propagate(fitted@control,length(range)^length(which))
          
           sq=list(range)
           sq=do.call('expand.grid',sq[rep(1,length(which))])
       
           for (i in seq(length(which))){
               fitted@control[which[i],'val']=     params(fitted)[which[i]]*sq[,i]
               fitted@control[which[i],'min']=min(fitted@control[which[i],'val'])*range[1]
               fitted@control[which[i],'max']=max(fitted@control[which[i],'val'])*range[2]}

           fitted@control[which,'phase']=-1
           }
        else
          fitted@control=profileGrid(fitted@control,which,range)
          
        if (!run) return(fitted)
        
        res=fit(fitted,index)
        
        rtn=fn(res)
        
        return(rtn)})

# # CIs
# cis <- max(surface) - qchisq(ci, 2)
# 
# do.call('contour', list(x=sort(profiled[[1]]), y=sort(profiled[[2]]), z=surface,
#                             levels=cis, add=TRUE, labcex=0.8, labels=ci))
# 
  
profileGrid=function(object,which,range=seq(0.95,1.05,length.out=11)){
    
    res=maply(range,function(x,which,ctl) {
      ctl[which,'val']=ctl[which,'val']*x
      ctl}, which=which, ctl=object)
    
    res=aperm(res,c(2:4,1))
    
    dmns=dimnames(object)
    dmns$iter=seq(dim(object)[3]*length(range))
    
    object=FLPar(array(c(res),unlist(lapply(dmns,length)),dmns))
    
    return(object)}
  
  