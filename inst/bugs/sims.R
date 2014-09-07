library(biodyn)
library(ggplotFL)
library(plyr)
library(diags)

load("/home/laurie/Desktop/sims.RData")

dat=subset(plot(sims)$data,qname!="Yield")[,c(1,3,7,8,9)]
dat=cast(dat,.id+year+qname~iter,value="data")


ggplot(subset(dat,.id!="Best Fit"))+
  geom_line(aes(year,`0.5`))+
  geom_ribbon(aes(year,ymin=`0.05`,ymax=`0.95`),alpha=.5)+
  geom_line(aes(year,`0.5`),subset(dat,.id!="Best Fit"),col="red")+
  facet_grid(qname~.id,scale="free_y")


load("/home/laurie/Desktop/bdHat.RData")

plot(bdHat)

harvest=FLQuant(c(fmsy(bdHat)),dimnames=list(year=20:50))
bdFwd  =fwd(bdHat,harvest=harvest)
plot(bdFwd)

catch  =FLQuant(c( refpts(bdHat)["msy"]),dimnames=list(year=20:50))
bdFwd  =fwd(bdHat,catch=catch)
plot(bdFwd)

stock  =FLQuant(c(bmsy(bdHat)),dimnames=list(year=20:50))
bdFwd  =fwd(bdHat,stock=stock)
plot(bdFwd)

bdFwds  =fwd(bdHat,FLQuants(harvest=harvest,catch=catch,stock=stock))
plot(bdFwds)

pe   =rlnorm(1,FLQuant(0,dimnames=list(year=1:50)),.2)
bdFwd=fwd(bdHat,harvest=harvest(bdHat)[,-dims(bdHat)$year],pe=pe)
plot(bdFwd)

harvest=rlnorm(100,log(harvest(bdHat))[,-dims(bdHat)$year],.1)
bdFwd=fwd(bdHat,harvest=harvest)
plot(bdFwd)

TAC=FLQuants(catch=catch*0.0,
             catch=catch*0.2,
             catch=catch*0.4,
             catch=catch*0.6,
             catch=catch*0.8,
             catch=catch*1.0)

bdTAC=fwd(bdHat,TAC)
names(bdTAC)=paste("MSY",c("0%","20%","40%","60%","80%","100%"))
plot(bdTAC)

