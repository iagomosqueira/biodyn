setGeneric("jackSummary", function(object, ...)
  standardGeneric("jackSummary"))

setMethod("jackSummary", signature(object="FLQuant"),
  function(object,sim,...) {
            
    n   <- dims(sim)$iter
            
    mn  <- object
    u   <- sim
    mnU <- apply(u, 1:5, mean)   
            
    SS  <- apply(sweep(u, 1:5, mnU,"-")^2, 1:5, sum)
            
    bias<- (n - 1) * (mnU - mn)
    se  <- sqrt(((n-1)/n)*SS)
            
   return(FLQuants(mean=mn, se=se, bias=bias))
   })

setMethod("jackSummary", signature(object="FLPar"),
  function(object,sim,...) {
            
   nms <-names(dimnames(object))
   idx <-seq(length(nms))[nms != 'iter']
   n   <-dims(sim)$iter 
            
   mn  <-object
   u   <-sim
   mnU <-apply(u, idx, mean)   
            
   SS <-apply(sweep(u, idx, mnU,"-")^2, idx, sum)
            
   bias <- (n - 1) * (mnU - mn)
   se   <- sqrt(((n-1)/n)*SS)
   
   cov  <-FLPar(cov(model.frame(u)[,dimnames(u)$params])*(n-1)*(n-1)/n)
   
   return(FLPars(hat=mn, mean=mnU, se=se, bias=bias, cov=cov))})
