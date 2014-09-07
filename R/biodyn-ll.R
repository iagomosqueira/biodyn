utils::globalVariables('calcSigmaFLQ')

# calcSigma
calcSigma=function(obs,hat=rep(0,length(obs)),error='log'){
  SS   <-sum((obs-hat)^2,na.rm=T)
  
  return((SS/length(hat))^.5)}

# calcLogLik
calcLogLik<-function(obs,hat=rep(0,length(obs)),error='log',type=1){
  
  logl<-function(se,obs,hat){
    SS<-sum((obs-hat)^2)
    
    n   <-length(obs)
    res <-(log(1/(2*pi))-n*log(se)-SS/(2*se^2))/2
    
    return(res)}
  
  se<-calcSigma(obs,hat,error=error)
  
  if (type==1) return(logl(se,obs,hat)) else
    if (type==2) return(-sum(dnorm(obs, hat, se, log=(error=='log')))) else
      if (type==3) return(sum((obs-hat)^2))}


llSigma=function(obs,hat=obs*0,dims=c(1,3:6)){
  
  hat=hat[,dimnames(obs)$year]
  
  SS =apply(obs-hat, dims, function(x) sum(x^2,na.rm=T))
  n  =apply(obs-hat, dims, function(x) sum(!is.na(x)))
  
  return((SS/n)^.5)}

llQ=function(obs,hat,dims=c(1,3:6),error='log'){
  
  yrs=dimnames(obs)$year[dimnames(obs)$year %in% dimnames(hat)$year]
  obs=obs[,yrs]
  
  res=switch(error,
             normal={q    =  apply(hat*obs,dims, function(x) sum(x))
                     q    =q/apply(hat,    dims, function(x) sum(hat*hat))
                     sigma=calcSigmaFLQ(obs/(q%*%hat))
                     
                     FLQuants(q=q,sigma=sigma)},
             log   ={q    =apply(log(obs)-log(hat), dims, 
                                 function(x) exp(sum(x,na.rm=T)/sum(!is.na(x))))
                     sigma=llSigma(log(obs),log(q%*%hat))
                     
                     FLQuants(q=q,sigma=sigma)},
             cv   ={res   =apply(obs/hat, dims, sum) #bug!
                    sigma2=llSigma(res)
                    q     =(-res+(res^2+4*apply(obs,dims,function(x) sum(!is.na()))*sigma2*
                                    apply(obs/hat, dims, function(x) sum((x)^2)))/
                              (2*apply(obs,dims,function(x) sum(!is.na()))*sigma2))          
                    
                    FLQuants(q=q,sigma=sigma)})
  
  return(res)}

# calcLogLik
calcLl<-function(obs,hat=obs*0,error='log',type=1){
  
  hat=hat[,dimnames(obs)$year]
  
  logl<-function(se,obs,hat=obs*0,dims=c(1,3:6)){
    
    SS  =apply(obs-hat, dims, function(x) sum((x)^2,na.rm=T))    
    n   =apply(obs, dims, function(x) sum(!is.na(x)))
    
    res =(log(1/(2*pi))-n*log(se)-SS/(2*se^2))/2
    
    return(res)}
  
  se=llSigma(obs,hat)
  
  if (type==1) return(logl(se,obs,hat))
  if (type==2) return(apply(obs-hat, dims, 
                      function(x) {
                           se=llSigma(x)  
                          -sum(dnorm(x, 0, se, log=(error=='log')))}))
  if (type==3) return(apply(obs-hat, dims, function(x) sum(x^2)))

  }

