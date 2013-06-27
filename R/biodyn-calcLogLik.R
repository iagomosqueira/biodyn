# calcSigma
calcSigma=function(obs,hat=rep(0,length(obs)),error="log"){
  SS   <-sum((obs-hat)^2,na.rm=T)
  
  return((SS/length(hat))^.5)}

# calcLogLik
calcLogLik<-function(obs,hat=rep(0,length(obs)),error="log",type=1){
  
  logl<-function(se,obs,hat){
    SS<-sum((obs-hat)^2)
    
    n   <-length(obs)
    res <-(log(1/(2*pi))-n*log(se)-SS/(2*se^2))/2
    
    return(res)}
  
  se<-calcSigma(obs,hat,error=error)
  
  if (type==1) return(logl(se,obs,hat)) else
    if (type==2) return(-sum(dnorm(obs, hat, se, log=(error=="log"), na.rm=TRUE))) else
      if (type==3) return(sum((obs-hat)^2))}
