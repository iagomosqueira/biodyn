quantiles<-function(x,probs=seq(0, 1, 0.25), 
                      na.rm= FALSE,
                      names= TRUE, 
                      type = 7){
    require(plyr)
  
    llply(res,quantile,probs=probs,na.rm=na.rm,names=names,type=type)}

FLQuantPoints<-function(x,...){
  require(plyr)
  
  res=FLQuants(x,...)
  
  res=llply(res,FLQuantPoint)
  return(res)}