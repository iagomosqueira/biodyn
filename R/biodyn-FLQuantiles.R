FLQuantiles<-function(x,...,
                      probs=seq(0, 1, 0.25), 
                      na.rm= FALSE,
                      names= TRUE, 
                      type = 7){
    require(plyr)
  
    res=FLQuants(x,...)

    res=llply(res,quantile,probs=probs,na.rm=na.rm,names=names,type=type)
    
    return(res)}

FLQuantPoints<-function(x,...){
  require(plyr)
  
  res=FLQuants(x,...)
  
  res=llply(res,FLQuantPoint)
  return(res)}