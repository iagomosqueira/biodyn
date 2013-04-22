#' biodyn Class
#' @description Creates an object of the biodyn class representing a biomass dynamic stock assessment model.
#' @name biodyn
#' @param model a factor or string that specifies the model type, has to be one of "fox", "schaefer", "pellat", "gulland", "fletcher", "shepherd", "logistic", "genfit"
#' @param params model parameters
#' @return biodyn object
#' @export
#' @examples bd=biodyn("pellat",FLPar(r=0.6,k=50000,p=1,b0=1))
setGeneric('biodyn',   function(model,params,...)  standardGeneric('biodyn'))
setMethod('biodyn', signature(model='factor',params="FLPar"),
          function(model,params,min=0.1,max=10,msy=NULL,catch=NULL,stock=NULL,...){
     
            args = list(...)
          
            dimnames(params)$params=tolower(dimnames(params)$params)

            if (!("b0" %in%  dimnames(params)$params)) 
              params=rbind(params,propagate(FLPar("b0"=1),dims(params)$iter))
            
            if (model=="pellat" & !("p" %in%  dimnames(params)$params)) 
              params=rbind(params,propagate(FLPar("p"=1),dims(params)$iter))
            
            if (!("k" %in%  dimnames(params)$params) & !is.null(msy))
              if (model=="pellat") params=rbind(params,"k"=FLPar(calcK(msy,params)))
          
            if (model=="pellat")
                params=params[c("r","k","p","b0"),]
            res        =biodyn()
            res@model  =model
            res@params =params 
            
            if (!is.null(stock))
               res@stock[]=params(res)["k"]*params(res)["b0"]
            
            if (!is.null(catch)){
               res@catch=catch
               res=fwd(res,catch=catch)
               }
            else  if (!is.null(stock)) {
               res@catch=window(res@stock,end=dims(res@catch)$maxyear-1)
               res@catch[]=NA}
            
            res@control=propagate(res@control,dims(params)$iter)
            nms=dimnames(res@control)$param[dimnames(res@control)$param %in% dimnames(res@params)$param]
            res@control[nms,  "val"]=res@params[nms,]
            res@control[nms,  "min"]=res@params[nms,]*min
            res@control[nms,  "max"]=res@params[nms,]*max
            
            if (!("b0" %in% nms))
               res@control["b0",c("min","max","val")]=c(0.75,1,1)
                
            # Load given slots
            for(i in names(args))
              slot(res, i) = args[[i]]
                          
            return(res)})

setMethod('biodyn', signature(model='character',params="FLPar"),
          function(model,params,min=0.1,max=10,msy=NULL,catch=NULL,stock=NULL,...) 
            biodyn(model=factor(model,levels=biodyn:::models),params,min=min,max=max,msy=msy,catch=catch,stock=stock,...))

setMethod('biodyn', signature(model='factor',params="missing"),
          function(model,params,min=min,max=max,msy=msy,catch=NULL,stock=NULL,...){
            
            args = list(...)
            
            res        =biodyn()
            res@model  =model
            
            nms=c(biodyn:::modelParams(model),"b0")
            par=rep(NA,length(nms))
            names(par)=nms
            
            res@params =FLPar(par) 
            res@params["b0"]=1
            
            if (!is.null(stock))
              res@stock[]=params(res)["k"]*params(res)["b0"]
            
            if (!is.null(catch))
              res@catch=catch
            else  if (!is.null(stock)) {
              res@catch=window(res@stock,end=dims(res@catch)$maxyear-1)
              res@catch[]=NA}
            
            # Load given slots
            for(i in names(args))
              slot(res, i) = args[[i]]
            
            return(res)})

setMethod('biodyn', signature(model='character',params="missing"),
          function(model=model,min=0.1,max=10.0,msy=NULL,catch=NULL,index=NULL,stock=NULL,...) 
            biodyn(model=factor(model,levels=biodyn:::models),min=min,max=max,catch=catch,stock=stock,...))

setMethod('biodyn', signature(model='missing',params="missing"),
          function(model,params,min=0.1,max=10.0,catch=catch,stock=stock,...) {
            args = list(...)
               
            res=new("biodyn")
            
            # Load given slots
            for(i in names(args))
              slot(res, i) = args[[i]]
            
          return(res)})

#' Checks class type
#'
#' @description Returns TRUE if object is of type biodyn
#' @param x biodyn class
#' @return TRUE or FALSE
#' @export
#' @examples
#' is.biodyn(biodyn()) 
is.biodyn = function(x)
  return(inherits(x, "biodyn"))

