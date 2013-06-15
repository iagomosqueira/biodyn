delta=function(par,object,fn){
    
    params(object)[dimnames(par)$params]=par
    object=fwd(object,catch=catch(object))
  
    rtn=fn(object)
    
    print(par)
    print(rtn)
    
    return(rtn)}
  
#   object=swon[[1]]
#   par=params(object)
#   
#   delta(par,object,fn=function(x) stock(x)[,ac(range(x)["maxyear"])]%/%bmsy(x))
# 
# object=fit(as(swon[[1]],"biodyn"),index(swon[[1]],F))
# hess  =hessian(func=delta,x=params(object)["r"],object=object, fn=function(x) stock(  x)[,ac(range(x)["maxyear"])]%/%bmsy(x))
# t.=data.frame(c(hess),c(vcov(object)))
# sum(t.[,1]*t.[,2])^.5
# 
# hess2 =hessian(func=delta, x=params(object),object=object, fn=function(x) harvest(x)[,ac(range(x)["maxyear"])]%/%fmsy(x))
# t.=data.frame(c(hess2),c(vcov(object)))
# sum(t.[,1]*t.[,2])^.5

