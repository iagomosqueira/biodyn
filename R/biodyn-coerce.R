#' as("biodyn", "FLStock")
#'
#' @description
#' Coerces an \code{FLStock} into a \code{biodyn} object.
#'
#' @name as
#' @family biodyn
#'
#' @rdname as
FLStock2biodyn=function(from){
  
  res      =biodyn()
  res@catch=catch(from)
  res@stock=window(stock(from),end=dims(from)$maxyear+1)
  
  dmns=dimnames(res@params)
  
  res@params=FLPar(as.numeric(NA),dimnames=dmns)
  res@params[]=c(.6,4*mean(res@catch),1,1)
  
  res@control[,'val']=res@params
  res@control[,'min']=res@control[,'val']*.1
  res@control[,'max']=res@control[,'val']*10
  res@control[,'phase']=c(1,1,-1,-1)
  
  range(res)[]=range(from)[c('minyear','maxyear')]
  
  res}

setAs('biodyn', 'FLStock',
      function(from) FLStock2biodyn(from))

FLBRP2biodyn=function(from){
  
  ## age based regference points
  msy =c(from@refpts["msy",   "yield"])
  bmsy=c(from@refpts["msy",   "biomass"])
  k   =c(from@refpts["virgin","biomass"])
  
  # bbiomass based reference points
  p   =optimise(function(p,bmsy,k) 
    (bmsy-k*(1/(1+p))^(1/p))^2, c(0.001,5), bmsy=bmsy,k=k)$minimum
  k=bmsy/((1/(1+p))^(1/p))
  r=msy/(k*(1/(1+p))^(1/p+1))
  b0=mean(stock(popl)[,1]/k)
  
  bd=biodyn()
  
  bd@params=FLPar(c(r=r,k=k,p=p,b0=b0))
  
  bd@catch=catch.obs(from)
  bd@stock=stock.obs(from)
  
  bd}

setAs('biodyn', 'FLBRP',
      function(from) FLBRP2biodyn(from,to))

setMethod('biodyn', signature(model='FLStock',params='missing'),
    function(model){
                  
      res      =new('biodyn')
      res@catch=catch(model)
      res@stock=window(stock(model),end=dims(model)$maxyear+1)
      
      dmns=dimnames(res@params)
      
      res@params=FLPar(as.numeric(NA),dimnames=dmns)
      res@params[]=c(.6,4*mean(res@catch),1,1)
      
      res@control[,'val']=res@params
      res@control[,'min']=res@control[,'val']*.1
      res@control[,'max']=res@control[,'val']*10
      res@control[,'phase']=c(1,1,-1,-1)
      
      range(res)[]=range(model)[c('minyear','maxyear')]
      
      res})

setGeneric('as.biodyn',   function(x,...)     standardGeneric('as.biodyn'))
setMethod('as.biodyn',signature(x='FLStock'),function(x){FLStock2biodyn(x)})
          