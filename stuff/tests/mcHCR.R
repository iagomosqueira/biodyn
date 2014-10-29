library(biodyn)
library(plyr)
library(ggplot2)

pe=rlnorm(500,FLQuant(0,dimnames=list(year=1:50)),0.5)

bd=window(sim(),end=30)
bd.=bd
bd@stock =propagate(bd@stock, 500)
bd=fwd(bd,harvest=harvest(bd)[,2:30],pe=pe)

for (i in seq(30,48,1))
  bd=fwd(bd,catch=hcr(bd,refYrs=i,yrs=i+1,tac=T,bndTac=c(0.9,1.1))$tac,pe=pe)

plot(bd)

trks=biodyn::kobe(bd,what="trks")
trks=mdply(data.frame(Year=seq(33,49,3)), function(Year) subset(trks,year<=Year))

pts =transform(biodyn::kobe(bd,what="pts",year=seq(33,49,3)),Year=year)[,c("stock","harvest","Year")]

kobePhase()+    
    geom_line(aes(stock,harvest),data=biodyn:::hcrPlot(bd.),col="brown",size=1.5)+    
    geom_path( aes(stock,harvest),data=subset(trks,pctl=="50%"))+
    geom_point(aes(stock,harvest),data=subset(pts,Year>=33))  +
    facet_wrap(~Year)
