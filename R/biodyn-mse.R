utils::globalVariables(c('br','srDev'))

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
mseBiodyn<-function(om,br,srDev,ctrl,prrs,
                    start,          end,         interval=3,
                    ftar  =0.70,    btrig =0.60,
                    fmin  =0.01,    blim  =0.01,
                    bndF  =NULL,
                    maxF  =2.0,     
                    uCV   =0.3,     phaseQ=1,   fishDepend=TRUE,
                    cmdOps=paste('-maxfn 500 -iprint 0 -est')){
 
  ## Get number of iterations in OM
  nits=c(om=dims(om)$iter, br=dims(params(br))$iter, rsdl=dims(srDev)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in OM")
  if (nits['om']==1) stock(om)=propagate(stock(om),max(nits))

  ## Cut in capacity
  maxF=mean(apply(fbar(window(om,end=start)),6,max)*maxF)
  
  ## hcr without feedback
  #  hcr parameters
#   par=hcrParam(btrig=btrig*br['msy','biomass'],
#                 blim =blim *br['msy','biomass'],
#                 ftar =ftar *br['msy','harvest'],
#                 fmin =fmin *br['msy','harvest'])
#   # loop over years
#   for (iYr in seq(start+rcvPeriod,range(om,'maxyear')-interval,interval)){
#     # data from last year then use to set TAC start of next year
#     TAC=tacFn(om,hcrFn(om,par,iYr-1,1+seq(interval),br,srDev)
#     om =fwd(om,catch=TAC,sr=br,sr.residuals=srDev)}
#   # save for reference  
  mou=om
   
  #### Observation Error (OEM) setup #######################
  ## Random variation for CPUE  
  cpue=oem(window(om,end=start),uCV)
  
  ## Loop round years
  mp =NULL
  hcr=NULL
  for (iYr in seq(start,range(om,'maxyear')-interval,interval)){
    #iYr = seq(start+rcvPeriod,range(om,'maxyear')-interval,interval)[1]
    cat('\n===================', iYr, '===================\n')
    
    ## use data from last year
    cpue=window(cpue,end=iYr-1)
    cpue[,ac(iYr-(interval:1))]=oem(om[,ac(iYr-(interval:1))],uCV)
    
    #### Management Procedure
    ## Set up assessment parameter options
    bd=biodyn:::biodyn(window(om,end=iYr-1))
    params(bd)[dimnames(ctrl)$param]=ctrl[dimnames(ctrl)$param,'val']
    
    bd@priors=prrs
    setParams( bd)=cpue
    setControl(bd)=params(bd)
    bd@control[dimnames(ctrl)$params,'phase'][]=ctrl[dimnames(ctrl)$params,'phase']
    bd@control['q1','phase']=phaseQ
    bd@control['q1','val']  =1

    ## fit
    bd =biodyn:::fit(bd,cpue,cmdOps=cmdOps)
    bd =biodyn:::fwd(bd,catch=catch(om)[,ac(iYr)])
        
    ## HCR
    hcrPar=hcrParam(ftar =ftar *fmsy(bd),
                     btrig=btrig*bmsy(bd),
                     fmin =fmin *fmsy(bd), 
                     blim =blim *bmsy(bd))
    hcrOutcome=biodyn:::hcr(bd,hcrPar,
                   hcrYrs=iYr+seq(interval),
                   bndF=bndF,
                   tac =TRUE)
            
    ## TACs for next year (iYtr+1) for n=interval years
    TAC  =hcrOutcome$tac
    TAC[]=rep(apply(TAC,6,mean)[drop=T],each=interval)
    
    #### Operating Model Projectionfor TAC
    om =fwd(om,catch=TAC,maxF=maxF,sr=br,sr.residuals=srDev)  

    #### Summary Statistics
    ## HCR actions, i.e. is biomass<Btrig?, what is F?, ..
    hcr =rbind(hcr,data.frame(yearHcr=min(as.numeric(dimnames(hcrOutcome$hvt)$year)),
                              #yearAss=rep(range(bd)[2],dims(bd)$iter),
                              model.frame(           hcrPar,drop=T)[,-5],
                              tac    =as.data.frame(apply(hcrOutcome$tac,6,mean),drop=T)[,'data'],
                              harvest=as.data.frame(apply(hcrOutcome$hvt,6,mean),drop=T)[,'data'],
                              stock  =as.data.frame(hcrOutcome$stock,drop=T)[,2]))
    
    ## Assessment parameters and reference points
    mp =rbind(mp,cbind(cbind(year=iYr,model.frame(params(bd))),
                       model.frame(refpts(bd))[,-4],
                       hcr))
    }
  
  ## save OM, projection without feedback, last assessment and MP summary
  return(list(om=om,mou=mou,bd=bd,mp=mp,oem=mcf(FLQuants(cpue=cpue,catch=catch(om)))))}


hcrFn=function(om,btrig,blim,ftar,fmin,start,end,interval,lag=seq(interval)){    
  
  a=(ftar-fmin)/(btrig-blim)
  b=ftar-a*btrig
  
  for (iYr in seq(start+rcvPeriod,range(om,'maxyear')-interval,interval)){
    stk=FLCore:::apply(stock(om)[,ac(iYr-1)],6,mean)
    
    trgt=(stk%*%a)+b  
    
    for (i in seq(dims(om)$iter)){
      FLCore:::iter(trgt,i)[]=max(FLCore:::iter(trgt,i),FLCore:::iter(fmin,i))
      FLCore:::iter(trgt,i)[]=min(FLCore:::iter(trgt,i),FLCore:::iter(ftar,i))} 
    
    dmns     =dimnames(trgt)
    dmns$year=as.character(iYr+lag)
    
    trgt=FLQuant(rep(c(trgt),each=length(lag)),dimnames=dmns)
    
    om=fwd(om,f=trgt,sr=br,sr.residuals=srDev)
  }
  
  return(om)}


