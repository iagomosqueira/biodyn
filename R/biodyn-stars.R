
ROregimeSHFT<-function(regLN,sig,series,shift=0){
  Shifts<-0
  startyr <-1
  if(shift!=0) {startyr<-shift}
  for (i in seq(startyr,(length(series)-(regLN)) ) )
  {
    window<-series[i:(regLN+(i-1))]
    basemean<-mean(window,na.rm=T)
    basevar<-var(window,na.rm=T)
    val<-sig*sqrt(2*basevar/regLN)
    
    # the below line would be changed to something like:
    # if (series[i+1] > basemean+val | series[i+1] <basemean-val)
    # for the original algorithm (I have not valed this)
    
    if (series[regLN+i] > basemean+val | series[regLN+i] <basemean-val)
    { 
      RSI<-0
      if (series[regLN+i] > basemean+val)
      {
        newmean<-basemean+val
        RSI<-(series[regLN+i]-newmean)/(regLN*(sqrt(basevar)))
        down<-0} else
        {
          newmean<-basemean-val
          RSI<-(newmean-series[regLN+i])/(regLN*(sqrt(basevar)))
          down<-1
        }
      
      origSIGN<-sign(RSI)
      counter<-min(regLN-1,(length(series)-(regLN+(i-1))))
      RSI<-0
      for (j in 1:counter)
      { 
        Shifts<-0
        if(down==0)
        {RSI<- RSI + (series[regLN+i+j-1]-newmean)/(regLN*(sqrt(basevar)))}else
        {RSI<- RSI + (newmean-series[regLN+i+j-1])/(regLN*(sqrt(basevar)))}
        
        
        if (sign(RSI)!=origSIGN | RSI == Inf) {break}
        if (j == counter)
        {
          Shifts<-i
          break
        }
      }
    }
    if(Shifts!=0) {break}
  }
  return(Shifts)
}


rod=function(val,year=NULL,ggplot=TRUE){
  #==set the assumed minimum regime length==
  rgLN<-10
  
  #== 'sig' should be based on an appropriate t-distribution given
  #== the number of observations in the series (but I didn't do that here)
  #== this significance level is roughly equal to p<0.1
  sig<-1.68
  
  #==iteration for finding all the shifts in a time series==
  Shift<-rep(NA,100)  						# vector for recording shifts
  if(length(val)>rgLN)
  {
    try(Shift[1]<-ROregimeSHFT(regLN=rgLN,sig=1.68,series=val),TRUE)
    counts<-2
    if(is.na(Shift[counts-1])==FALSE) 
    {
      while(Shift[counts-1]>0 & (Shift[counts-1]+rgLN+rgLN-1)<length(val))
      {
        Shift[counts]<-ROregimeSHFT(rgLN,sig,val,(Shift[counts-1]+rgLN-1))
        counts<-counts+1
      }
    }
  }   
  
  #===plot results=== 
  
  ShiftsVec<-Shift[!is.na(Shift)]
  
  #plot(val,type='l')
  mn=NULL
  sd=NULL
  ln=NULL
  
  for(i in 1:length(ShiftsVec))
  {
    if(i==1){  
      regMean<-mean(val[1:(Shift[i]+rgLN-1)])
      regSD<-sd(val[1:(Shift[i]+rgLN-1)])
      mn=c(mn,regMean)
      sd=c(sd,regSD)
      ln=c(ln,Shift[i]+rgLN-1)
      
#       polygon(x=c(1,Shift[i]+rgLN-1,Shift[i]+rgLN-1,1),
#               y=c(regMean+regSD,regMean+regSD,regMean-regSD,regMean-regSD),
#               border=NA,col='#0000ff55')
    }
    
    if(i>1){
      endInd<-ShiftsVec[i]+rgLN-1
      if(is.na(ShiftsVec[i+1]))
      {endInd<-length(val)}
      
      regMean<-mean(val[(Shift[i-1]+rgLN):endInd])
      regSD  <-sd(  val[(Shift[i-1]+rgLN):endInd])
      mn=c(mn,regMean)
      sd=c(sd,regSD)
      ln=c(ln,endInd)
      
#       polygon(x=c(Shift[i-1]+rgLN,
#                   endInd,
#                   endInd,         
#                   Shift[i-1]+rgLN),
#               y=c(regMean+regSD,
#                   regMean+regSD,
#                   regMean-regSD,
#                   regMean-regSD),
#               border=NA,col='#0000ff55')
    }
  }
  
  #print(ln)
  if (!is.null(year)){
    minyear=c(year[1],year[rev(rev(ln)[-1])]+1)
    maxyear=min(year)+ln-1
    
    rtn=data.frame(mn=mn,sd=sd,i=factor(seq(length(mn))),
                   ln=ln,minyear=minyear,maxyear=maxyear)
    
    if (!ggplot) return(rtn)
    
    forGG=function(dat) data.frame(i=dat$i,
                                   x=with(dat,c(minyear,minyear,maxyear,maxyear)),
                                   y=with(dat,c(mn+sd,  mn-sd,  mn-sd,  mn+sd)))
    
    return(forGG(rtn))
  }
  else
    return(data.frame(mn=mn,sd=sd,i=factor(seq(length(mn))),
                      ln=ln))

}