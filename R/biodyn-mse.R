  mseBiodyn<-function(OM,start,
                      brp,srRsdl=FLQuant(1,dimnames=dimnames(window(rec(OM),start=start))),
                      CV   =0.1,
                      Ftar =0.75,    Btrig=0.75,
                      Fmin =Ftar*0.1,Blim =Btrig*0.0001,
                      Bpct =0.5,     Fpct =0.5,
                      jk    =FALSE,
                      bounds=NULL,
                      fishDepend=TRUE){
  
  ## Get number of iterations in OM
  nits=c(OM=dims(OM)$iter, sr=dims(params(brp))$iter, rsdl=dims(srRsdl)$iter)
  if (length(unique(nits))>=2 & !(1 %in% nits)) ("Stop, iters not '1 or n' in OM")
  nits=max(nits)
  stock(OM)=propagate(stock(OM),nits)

  #### Observation Error (OEM) setup #######################
  ## Random variation for Catch & CPUE, CV=0.25%
  bd       =biodyn(OM)
  params(bd)["k"]= refpts(brp)["virgin","biomass"]*2
 
  params(bd)["k"]=5.0e+09
    
  params(bd)["r"]= log(lambda(leslie(brp,c(refpts(brp)["crash","harvest"]))))
  params(bd)["b0"]=.9

  params(bd)["r"]=0.3
    
  bd@mng=FLPar(a=1) 
  bd       =propagate(bd,nits)

  if (fishDepend) cpue=catch(OM)/fbar(OM) else cpue=stock(OM)
  cpue=cpue/mean(cpue,na.rm=T)
  cpue=cpue/sd(cpue,na.rm=T)
  cpue=cpue*rlnorm(dim(cpue)[6],0,CV)
  setParams( bd)<-cpue
  setControl(bd)=params(bd)
  control(bd)[c("p","b0"),"phase"][]=-1
  
  if (!is.null(bounds)) bd@bounds=bounds
  ## Loop round years
  
  for (iYr in start:(range(OM,"maxyear")-2)){
  #iYr = (start:(range(OM,"maxyear")-2))[1]
     cat("===================", iYr, "===================\n")
     ## add year to FLBD
     bd                 =window(bd,  end=iYr)
     cpue               =window(cpue,end=iYr)
     catch(bd)[,ac(iYr)]=computeCatch(OM)[,ac(iYr)]
     
     if (jk){
      hv =hcrJK(bd,Ftar,Btrig,Fmin,Blim,Fpct,Bpct) 
     }else{ 
      bd =fit(bd,cpue)
 
      hv =hcr(bd,FLPar(Ftar=Ftar,Btrig=Btrig,Fmin=Fmin,Blim=Blim)) 
#     return(list(bd,cpue))
      }
     
     TAC=tac(bd,hv)
   
     #### Now you have TAC take it from OM
     ctrl    <-fwdControl(data.frame(year=iYr+2,max=c(NA,2),quantity=c("catch","f")))
     dms     <-dimnames(ctrl@trgtArray)
     dms$iter<-1:nits
     ctrl@trgtArray<-array(NA,lapply(dms,length),dms)
     ctrl@trgtArray[1,"val", ]<-TAC
     ctrl@trgtArray[2,"max", ]<-.99
   
     OM <-fwd(OM,ctrl=ctrl,sr=brp,sr.residuals=srRsdl)
     print(iYr)
     }

   return(list(OM=OM,MP=bd,cpue=cpue))}
