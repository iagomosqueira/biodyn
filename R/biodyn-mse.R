
mseBiodyn<-function(om,brp,ctrl,start=60,end=100,interval=3,
                    ftar =0.70,  btrig=0.60,
                    fmin =0.01,  blim =0.01,
                    seed=7890,   nits=100,
                    uCV   =0.1,  recCV=0.3,
                    fishDepend=TRUE,
                    bndF=c(1,Inf))
{
  
  set.seed(seed)
  om   =fwdWindow(om,end=end+interval,brp)
  srDev=rlnorm(nits,FLQuant(0,dimnames=dimnames(window(rec(om),end=end+interval))),recCV)
  
  om=fwd(om,f=fbar(om)[,2:(start+10)],sr=brp,sr.residuals=srDev)
  om=fwd(om,f=FLQuant(c(refpts(brp)["msy","harvest"]),dimnames=dimnames(fbar(om)[,(start+1):(end+interval)])),sr=brp,sr.residuals=srDev)
  
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, sr=dims(params(brp))$iter, rsdl=dims(srDev)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in OM")
  nits=max(nits)
  stock(om)=propagate(stock(om),nits)
  
  #### Observation Error (OEM) setup #######################
  ## Random variation for CPUE
  
  cpue=oem(window(om,end=start),uCV)
  
  ## cut in capacity
  maxF=max(fbar(window(om,end=start)))*.5
  
  ## Loop round years
  interval=3
  for (iYr in seq(start,range(om,"maxyear")-interval,interval)){
    #iYr = (start:(range(om,"maxyear")-2))[1]
    cat("===================", iYr, "===================\n")
    
    cpue=window(cpue,end=iYr)
    cpue[,ac(iYr-(interval:1)+1)]=oem(om[,ac(iYr-(interval:1)+1)],uCV)
    
    bd=biodyn(om)
    params(bd)[dimnames(ctrl)$param]=ctrl[dimnames(ctrl)$param,"val"]
    
    setParams( bd)=cpue
    setControl(bd)=params(bd)
    bd@control[c("p","b0","q1"),"phase"][]=-1
    bd@control["q1","val"]=1
    
    bd =fit(bd,cpue)
    
    hv =hcr(bd,hcrParams(ftar =ftar *fmsy(bd),
                         btrig=btrig*bmsy(bd),
                         fmin =fmin *fmsy(bd), 
                         blim =blim *bmsy(bd)),
            bndF=bndF)
    
    TAC=window(tac(bd,hv),end=iYr+interval)
    TAC[]=TAC[,1]    
    
    om <-fwd(om,catch=TAC,sr=brp,maxF=maxF,sr.residuals=srDev)
  }
  
  #plot(window(om,end=100))
  
  return(list(om=om,MP=bd,cpue=cpue))}


oem=function(om,cv,fishDepend=FALSE){
  
  nits=dims(catch(om))$iter
  rnd=rlnorm(nits,FLQuant(0,dimnames=list(year=dims(om)$minyear:dims(om)$maxyear)),cv)
  
  if (fishDepend) 
    cpue=rnd*catch(om)/fbar(om)
  else 
    cpue=rnd*computeStock(om)
  
  cpue}
