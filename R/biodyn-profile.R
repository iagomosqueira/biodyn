  ##############################################################
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
#' /dontrun{
#' library(aspic)
#' library(biodyn)
#' data(asp)
#' cpue=FLQuant(asp@index$index,dimnames=list(year=asp@index$year))
#' bd=as(asp,"biodyn")
#' setParams(bd)  =cpue
#' setControl(bd,.01,100) =params(bd)
#' res=profile(bd,cpue,which="r",fixed=c("b0","p"),range=c(.95,1.1),maxstep=51)
#' ggplot(res)+geom_line(aes(r,ll))
#' }

# ### debugging stuff
# data(bd)
# fitted=biodyn(factor("pellat"),params(bd),catch=catch(bd))
# cpue=rlnorm(1,log(stock(bd)),.2)[,-60]
# setParams(fitted)     =cpue
# 
# 
# attach(list(maxsteps=11, range=0.5, ci=c(0.25, 0.5, 0.75, 0.95),
#             plot=TRUE,fixed=c()))
# which="r"
# fixed=c("p","b0")
# ###
# res=profile(bd,which="r",fixed=c("b0","p"),cpue,range=c(1.2,3.0))
# ggplot(res)+geom_line(aes(r,ll))
# v <- ggplot(res, aes(r, k, z = ll))
# v <- ggplot(res, aes(r, k, z = ll))+ stat_contour(aes(colour = ..level..), size = 1)
setMethod("profile", signature(fitted="biodyn"),
      function(fitted,cpue,which,fixed=c(),
                   maxsteps=11, range=0.5,
                   fn   =function(x) cbind(model.frame(params(x)),model.frame(x@objFn)[,-3],model.frame(refpts(x))[,-4]),
                   run  =TRUE,...){
  
        if (dims(cpue)$maxyear>=dims(stock(fitted))$maxyear) stop("cpue years greater in length than stock")
        if (dims(fitted)$iter>1) stop("can only be done for a single iter")
                 
        #setControl(fitted)=params(fitted)
        control(fitted)=propagate(control(fitted),maxsteps^length(which))
          
        if (length(range)==1) range=c(range,2-range)
        
        sq=list(seq(range[1],range[2],length.out=maxsteps))
        sq=do.call("expand.grid",sq[rep(1,length(which))])
       
        for (i in seq(length(which))){
            control(fitted)[which[i],"val"]=     params(fitted)[which[i]]*sq[,i]
            control(fitted)[which[i],"min"]=min(control(fitted)[which[i],"val"])*range[1]
            control(fitted)[which[i],"max"]=max(control(fitted)[which[i],"val"])*range[2]}

        control(fitted)[c(fixed,which),"phase"]=-1
        
        if (!run) return(fitted)
        
        res=fit(fitted,cpue)
        
        rtn=fn(res)
        
        return(rtn)})

# # CIs
# cis <- max(surface) - qchisq(ci, 2)
# 
# do.call('contour', list(x=sort(profiled[[1]]), y=sort(profiled[[2]]), z=surface,
#                             levels=cis, add=TRUE, labcex=0.8, labels=ci))
# 
  