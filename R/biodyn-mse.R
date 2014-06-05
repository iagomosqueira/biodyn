## Observation Error Model
oem=function(om,cv,fishDepend=FALSE){
  
  nits=max(dims(stock(om))$iter,dims(catch(om))$iter)
  rnd=rlnorm(nits,FLQuant(0,dimnames=list(year=dims(om)$minyear:dims(om)$maxyear)),cv)
  
  if (fishDepend) 
    cpue=rnd*catch(om)/fbar(om)
  else 
    cpue=rnd*computeStock(om)
  
  cpue}

## MSE function
mseBiodyn<-function(om,brp,srDev,ctrl,
                    start,          end,         interval=3,
                    ftar  =0.70,    btrig =0.60,
                    fmin  =0.01,    blim  =0.01,
                    bndF  =NULL,
                    maxF  =2.0,     
                    uCV   =0.1,     phaseQ=-1,   fishDepend=TRUE,
                    seed =7890){
 
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, br=dims(params(brp))$iter, rsdl=dims(srDev)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in OM")
  if (nits["om"]==1) stock(om)=propagate(stock(om),max(nits))
  
  set.seed(seed)
  
  ## Cut in capacity
  maxF=mean(apply(fbar(window(om,end=start)),6,max)*maxF)
  
  ## hcr without feedback
  #  hcr parameters
#   par=hcrParams(btrig=btrig*brp["msy","biomass"],
#                 blim =blim *brp["msy","biomass"],
#                 ftar =ftar *brp["msy","harvest"],
#                 fmin =fmin *brp["msy","harvest"])
#   # loop over years
#   for (iYr in seq(start+rcvPeriod,range(om,"maxyear")-interval,interval)){
#     # data from last year then use to set TAC start of next year
#     TAC=tacFn(om,hcrFn(om,par,iYr-1,1+seq(interval),brp,srDev)
#     om =fwd(om,catch=TAC,sr=brp,sr.residuals=srDev)}
#   # save for reference  
  prj=om
   
  #### Observation Error (OEM) setup #######################
  ## Random variation for CPUE  
  cpue=oem(window(om,end=start),uCV)
  
  ## Loop round years
  mp =NULL
  hcr=NULL
  for (iYr in seq(start,range(om,"maxyear")-interval,interval)){
    #iYr = seq(start+rcvPeriod,range(om,"maxyear")-interval,interval)[1]
    cat("===================", iYr, "===================\n")
    
    ## use data from last year
    cpue=window(cpue,end=iYr-1)
    cpue[,ac(iYr-(interval:1))]=oem(om[,ac(iYr-(interval:1))],uCV)
    
    ## assessment 
    bd=biodyn(window(om,end=iYr-1))
    params(bd)[dimnames(ctrl)$param]=ctrl[dimnames(ctrl)$param,"val"]
    
    setParams( bd)=cpue
    setControl(bd)=params(bd)
    bd@control[dimnames(ctrl)$params,"phase"][]=ctrl[dimnames(ctrl)$params,"phase"]
    bd@control["q1","phase"]=phaseQ
    bd@control["q1","val"]  =1
    
    bd =fit(bd,cpue)
    bd =fwd(bd,catch=catch(om)[,ac(iYr)])
    
    ## HCR
    hcrPar=hcrParams(ftar =ftar *fmsy(bd),
                     btrig=btrig*bmsy(bd),
                     fmin =fmin *fmsy(bd), 
                     blim =blim *bmsy(bd))
    hcrOutcome=hcr(bd,hcrPar,
                   hcrYrs=iYr+seq(interval),
                   bndF=bndF,
                   tac =TRUE)
            
    ## Set TACs for next year (iYtr+1) for n=interval years
    TAC  =hcrOutcome$tac
    TAC[]=rep(apply(TAC,6,mean)[drop=T],each=interval)
    
    hcr =rbind(hcr,data.frame(yearHcr=rep(as.numeric(dimnames(hcrOutcome$hvt)$year),dims(bd)$iter),
                              #yearAss=rep(range(bd)[2],dims(bd)$iter),
                              model.frame(           hcrPar,drop=T),
                              tac    =as.data.frame(hcrOutcome$tac,drop=T)[,2],
                              harvest=as.data.frame(hcrOutcome$hvt,drop=T)[,2],
                              stock  =as.data.frame(hcrOutcome$stock,drop=T)[,2])[,-7])
    
    ## save assessment parameters and reference points
    mp =rbind(mp,cbind(cbind(year=iYr,model.frame(params(bd))),
                       model.frame(refpts(bd))[,-4],
                       hcr))
    
    om =fwd(om,catch=TAC,maxF=maxF,sr=brp,sr.residuals=srDev)  
    }
  
  ## save OM, projection without feedback, last assessment, assessed parameters and cpue
  return(list(om=om,prj=prj,bd=bd,mp=mp))}


hcrFn=function(om,btrig,blim,ftar,fmin,start,end,interval,lag=seq(interval)){    
  
  a=(ftar-fmin)/(btrig-blim)
  b=ftar-a*btrig
  
  for (iYr in seq(start+rcvPeriod,range(om,"maxyear")-interval,interval)){
    stk=FLCore:::apply(stock(om)[,ac(iYr-1)],6,mean)
    
    trgt=(stk%*%a)+b  
    
    for (i in seq(dims(om)$iter)){
      FLCore:::iter(trgt,i)[]=max(FLCore:::iter(trgt,i),FLCore:::iter(fmin,i))
      FLCore:::iter(trgt,i)[]=min(FLCore:::iter(trgt,i),FLCore:::iter(ftar,i))} 
    
    dmns     =dimnames(trgt)
    dmns$year=as.character(iYr+lag)
    
    trgt=FLQuant(rep(c(trgt),each=length(lag)),dimnames=dmns)
    
    om=fwd(om,f=trgt,sr=brp,sr.residuals=srDev)
  }
  
  return(om)}


