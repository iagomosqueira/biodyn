library(biodyn)
library(ggplotFL)
library(plyr)

bd=simBiodyn()

Us  =FLQuants("Unbiased"     =rlnorm(1,log((stock(bd)[,-dims(bd)$year]+
                                              stock(bd)[,-1])/2),0.2),
              "Increase in q"=rlnorm(1,log((stock(bd)[,-dims(bd)$year]+
                                              stock(bd)[,-1])/2),0.2))

plot(Us)

bd1=bd
bds=bd

setParams( bd1)=Us[[1]]
setControl(bd1)=params(bd1)

bd1@control[3:4,"phase"]=-1
bd1=fit(bd1,index=Us[[1]])
bd1@control[,c("min")]=bd1@params*0.1
bd1@control[,c("val")]=bd1@params
bd1@control[,c("max")]=bd1@params*10

prfl=profile(bd1,which='r',index=Us[[1]],
             range=seq(0.90,1.05,.001),fn=fn)
ggplot(subset(melt(prfl[,c("r","ll.u1")],id="r"),value<0))+
  geom_line(aes(r,value,group=variable,col=variable))+
  facet_wrap(~variable,scales="free")


setParams( bds)=Us
setControl(bds)=params(bds)

bds@control[3:4,"phase"]=-1
bds=fit(bds,index=Us)
bds@control[,c("min")]=bds@params*0.1
bds@control[,c("val")]=bds@params
bds@control[,c("max")]=bds@params*10

fn=function(x) cbind(model.frame(params(x)["r"]),
                     ll=model.frame(x@ll)[,-3],
                     ll=apply(x@ll,2,sum))

prfl=profile(bds,which='r',index=Us,
             range=seq(0.70,1.05,.001),fn=fn)
ggplot(subset(melt(prfl[,c("r","ll.u1","ll.u2","1")],id="r"),value<10))+
  geom_line(aes(r,value,group=variable,col=variable))+
  facet_wrap(~variable,scales="free",ncol=1)
