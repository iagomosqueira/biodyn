
## ----knitr_init, echo=FALSE, results="asis"------------------------------
library(knitr)

# output:
#   rmdformats::html_clean:
#     fig_width: 6
#     fig_height: 6
#     highlight: pygments

## Global options
options(max.print="75")
opts_chunk$set(fig.path="out/",
               echo =TRUE,
               cache=TRUE,
               cache.path="cache/",
               prompt=FALSE,
               tidy=TRUE,
               comment=NA,
               message=FALSE,
               warning=FALSE)

opts_knit$set(width=75)


## ------------------------------------------------------------------------
library(biodyn)


## ----,echo=FALSE---------------------------------------------------------
library(ggplotFL)
library(plyr)
library(diags)


## ----,eval=FALSE---------------------------------------------------------
## bd =biodyn()


## ----,eval=FALSE---------------------------------------------------------
## bd=biodyn(catch=FLQuant(100,dimnames=list(year=1990:2010)))


## ----,eval=FALSE---------------------------------------------------------
## data(ple4)
## bd =as(ple4,"biodyn")


## ----,eval=FALSE---------------------------------------------------------
## library(aspic)
## asp=aspic("http://http://rscloud.iccat.int/kobe/swon/2013/aspic/run2/")
## bd =as(asp,"biodyn")


## ------------------------------------------------------------------------
bd=simBiodyn()


## ----, fig.margin=TRUE, figure.width=4, figure.height=6, fig.cap="Simulated CPUE series"----
bd=simBiodyn()
bd=window(bd,end=49)


plot(bd)


## ----, fig.margin=TRUE, echo=FALSE, eval=FALSE, figure.height=4, fig.cap="Simulated CPUE series"----
## x=fwd(x,harvest=rlnorm(200,log(harvest(x)[,-1]),.2))
## plot(x,worm=3)


## ----, fig.margin=TRUE, echo=FALSE, eval=FALSE, fig.cap="Simulated CPUE series"----
## plotEql(x)


## ----, fig.margin=TRUE, fig.cap="Simulated CPUE series"------------------
library(reshape)
x=simBiodyn()
plotPrd(x)+
  geom_path( aes(stock,catch),model.frame(FLQuants(x,"stock","catch")))+
  geom_point(aes(stock,catch),model.frame(FLQuants(x,"stock","catch")))


## ----, fig.margin=TRUE, eval=FALSE, fig.cap="Simulated CPUE series"------
## plotMSE()


## ----, echo=TRUE, fig.margin=TRUE, fig.height=4, fig.cap="Simulated stock"----
bd=simBiodyn()


## ----, fig.margin=TRUE, fig.cap="Simulated CPUE series"------------------
cpue=(stock(bd)[,-dims(bd)$year]+
      stock(bd)[,-1])/2
cpue=rlnorm(1,log(cpue),.2)

ggplot(as.data.frame(cpue))+
  geom_point(aes(year,data))+
  geom_line(aes(year,data),data=as.data.frame(stock(bd)),col="salmon")


## ----,eval=FALSE---------------------------------------------------------
## bd=biodyn(catch=catch(bd),msy=mean(catch(bd)))


## ----,fig.margin=TRUE----------------------------------------------------
setParams( bd)=cpue
params(bd)


## ----,fig.margin=TRUE----------------------------------------------------
setControl(bd)=params(bd)
  
bd@control


## ----, eval=TRUE, fig.margin=TRUE, fig.height=6,fig.cap="A comparison of the true and fitted time series"----
bd@control[3:4,"phase"]=-1
bdHat=fit(bd,cpue)

# plot(biodyns("True"=bd,"Hat"=bdHat))+
#   theme(legend.position="bottom")
save(bdHat,cpue,file="/home/laurie/Desktop/bdHat.RData")


## ------------------------------------------------------------------------
params(bdHat)
params(bdHat)/params(bd)


## ----,echo=TRUE----------------------------------------------------------
rsdl=bdHat@diags

head(rsdl)


## -----rsdl5, fig.width=4,fig.height=4,fig.margin=TRUE, fig.cap="Quantile-quantile plot to compare residual distribution with the normal distribution."----
ggplot(rsdl)                                           +
  geom_point( aes(qqx,qqy))                            +
  stat_smooth(aes(qqx,qqHat),method="lm",se=T,fill="blue", alpha=0.1)         +
  theme_ms(14,legend.position="bottom")               


## -----rsdl2, fig.margin=TRUE, fig.height=4, figwidth=4, fig.cap="Observed CPUE verses fitted, blue line is a linear resgression fitted to points, black the y=x line."----
ggplot(with(rsdl, data.frame(obs=stdz(index),hat=stdz(hat))))+
          geom_abline(aes(0,1))                         +
          geom_point( aes(obs,hat))                     +
          stat_smooth(aes(obs,hat),method="lm", se=F)    +
          theme_ms(14,legend.position="bottom")            +
          xlab("Fitted") + ylab("Observed")


## -----rsdl3,fig.height=2, fig.margin=TRUE, fig.cap="Residuals by year, with lowess smoother"----
dat=transform(subset(rsdl,!is.na(residual), residual=stdz(residual,na.rm=T)))
ggplot(aes(year,residual),data=dat)  +
  geom_hline(aes(yintercept=0))      +
  geom_point()                       +
  stat_smooth(method="loess",se=F)  +
  theme_ms(14,legend.position="bottom") 


## -----rsdl6,fig.height=2, fig.margin=TRUE, fig.cap="Plot of residuals against fitted value, to check variance relationship."----
ggplot(aes(hat, residual),data=subset(rsdl,!is.na(hat) & !is.na(residual)))   +
  geom_hline(aes(yintercept=0))         +
  geom_point()                          +
  stat_smooth(method="loess",se=F)   +
  theme_ms(14,legend.position="bottom")


## -----rsdl4, fig.width=4,fig.width=4,fig.margin=TRUE, fig.cap="Plot of autocorrelation, i.e. $residual_{t+1}$ verses $residual_{t}$."----
ggplot(rsdl)                                              +
  geom_point( aes(residual,residualLag))                  +
  stat_smooth(aes(residual,residualLag),method="lm",se=F) +
  geom_hline(aes(yintercept=0))                           +
  xlab(expression(Residual[t])) + 
  ylab(expression(Residual[t+1])) +
  theme_ms(14,legend.position="bottom")  


## ----,eval=FALSE---------------------------------------------------------
## res=profile(bdHat,which='r',fixed=c('b0','p'),cpue,range=seq(0.97,1.03,.002))
## ggplot(res)+geom_line(aes(r,ll))


## ----,eval=FALSE---------------------------------------------------------
## res=profile(bdHat,which=c('r','k'),fixed=c('b0','p'),cpue,range=seq(0.97,1.03,.002))
## ggplot(res, aes(r, k, z=ll))+
##   stat_contour(aes(colour = ..level..), size = 1)


## ----,eval=FALSE---------------------------------------------------------
## library(biodyn)
## library(ggplot2)
## library(plyr)
## library(reshape)
## 
## bd=simBiodyn()
## bd=window(bd,end=49)
## cpue =(stock(bd)[,-dims(bd)$year]+
##       stock(bd)[,-1])/2
## cpue1=rlnorm(1,log(cpue),.2)
## cpue2=rlnorm(1,log(cpue),.2)
## setParams(bd) =FLQuants("1"=cpue1,"2"=cpue2)
## setControl(bd)=params(bd)
## 
## bd@control[3:4,"phase"]=-1
## 
## bd=fit(bd,index=FLQuants("1"=cpue,"2"=cpue2))
## 
## prfl=profile(bd,which='r',fixed=c('b0','p'),index=FLQuants("1"=cpue,"2"=cpue2),range=seq(0.97,1.03,.002))
## ggplot(res)+geom_line(aes(r,ll))


## ----,eval=FALSE,echo=FALSE----------------------------------------------
## bd2=fit(bdHat,cpue,cmdOps="-lprof")
## save(bd2,bdHat,cpue,file="/home/laurie/Desktop/bd2.RData")
## prf=subset(bd@profile, param %in% c("bbmsy","ffmsy"))
## prf=data.frame(What="Profile",t(daply(prf, .(param), with, sample(value,500,prob=p,replace=T))))
## names(prf)[2:3]=c("Stock","Harvest")


## ----, fig.margin=TRUE---------------------------------------------------
sims=biodyns("Best Fit"=bd)


## ----,figure.height=6, fig.margin=TRUE-----------------------------------
save(bdHat,cpue,file="/home/laurie/Desktop/bdHat.RData")
sims[["Vcov"]]=mvn(bdHat,500,nms=c("r","k"),fwd=TRUE)

plot(sims[["Vcov"]])

save(sims, file="/home/laurie/Desktop/sims.RData")


## ----, eval=FALSE, fig.height=4------------------------------------------
## cpueSim=bdHat@diags[,c("year","hat")]
## names(cpueSim)[2]="data"
## cpueSim=as.FLQuant(cpueSim)
## 
## #cv(diags["residuals"])
## 
## cpueSim=rlnorm(100,log(cpueSim),0.25)
## 
## cpueSim[cpueSim==0]=NA
## 
## plot(sims[["CPUE"]])
## 
## sims[["CPUE"]]=fit(propagate(bdHat,100),cpueSim)


## ----,fig.height=4,fig.margin=TRUE, fig.cap="Plot predicted stock trend by index"----
sims[["Jack Knife"]]=fit(bdHat,jackknife(cpue))


## ----,fig.height=4,fig.margin=TRUE, fig.cap="Plot predicted stock trend by index"----
plotJack(sims[["Jack Knife"]],bdHat)


## ----, fig.height=4, fig.margin=TRUE, fig.cap="Plot predicted stock trend by index"----
sims[["MCMC"]]=fit(bdHat,cpue,cmdOps=c("-mcmc 100000, -mcsave 50000"))

plot(sims[["MCMC"]])


## ----,fig.height=4,fig.margin=TRUE,--------------------------------------
acf(c(params(sims[["MCMC"]])["r"]))


## ----,fig.margin=TRUE----------------------------------------------------


## ------------------------------------------------------------------------
head(bdHat@mng)


## ----,eval=FALSE---------------------------------------------------------
## bdHat@mngVcov


## ----,fig.margin=TRUE----------------------------------------------------
currentState   =bdHat@mng[c("bbmsy","ffmsy"),"hat",drop=T]
currentStateVar=bdHat@mngVcov[c("bbmsy","ffmsy"),
                           c("bbmsy","ffmsy"),drop=T]

mvrnorm(10,currentState,currentStateVar)

ggplot(data=as.data.frame(mvrnorm(100,currentState,currentStateVar)),
       aes(bbmsy,ffmsy))+geom_point()+
  geom_hline(aes(yintercept=1))+
  geom_vline(aes(xintercept=1))+
  xlab(expression(B:B[MSY]))+
  ylab(expression(F:F[MSY]))


## ----,fig.margin=TRUE----------------------------------------------------
library(kobe)
 setGeneric('kobe',          function(object,method,...)    standardGeneric('kobe'))

source('~/Desktop/flr/pkgs/biodyn/R/biodyn-kobe.R')


## -----7------------------------------------------------------------------
df=kobe(sims)
ggplot(subset(df,year==49)) + 
  geom_density(aes(x=stock, y=..count..), position = "stack",fill="red") +
  geom_vline(aes(xintercept=1))          +
  facet_wrap(~.id,scales="free_y",ncol=1)
  #scale_x_continuous(limits=c(0.5,1.5))


## -----8------------------------------------------------------------------
ggplot(subset(df,year==49)) + 
  geom_density(aes(x=harvest, y=..count..), position = "stack",fill="red") +
  geom_vline(aes(xintercept=1))          +
  facet_wrap(~.id,scales="free_y",ncol=1)+
  scale_x_continuous(limits=c(0,1.5))


## -----9------------------------------------------------------------------
library(kobe)
kobePhase()+
  geom_point(aes(stock,harvest),data=subset(df,year==49))+
  facet_wrap(~.id,ncol=2)


## ----, fig.margin=TRUE, fig.cap=""---------------------------------------
harvest=rlnorm(1,log(harvest(bdHat))[,-dims(bdHat)$year],.1)
bdHat     =fwd(bdHat,harvest=harvest)

plot(bdHat)


## -----9b-----------------------------------------------------------------
bd   =simBiodyn()


## ----,fig.margin=TRUE----------------------------------------------------
## simulate HCRs, annual, tri-annula, F bound, TAC bound
bd=window(bd,end=29)
for (i in seq(29,49,1))
  bd=fwd(bd,harvest=hcr(bd,refYrs=i,yrs=i+1)$hvt)
simHCR=biodyns("1"=bd)

plot(bd)


## ----,fig.margin=TRUE----------------------------------------------------
bd=window(bd,end=29)
for (i in seq(29,49,3))
  bd=fwd(bd,harvest=hcr(bd,refYrs=i,yrs=i+1:3)$hvt)
simHCR[["3"]]=bd
save(simHCR,file="/home/laurie/Desktop/sims.RData")
plot(simHCR)


## ----,fig.margin=TRUE----------------------------------------------------
bd=window(bd,end=29)
for (i in seq(29,49,1))
  bd=fwd(bd,harvest=hcr(bd,refYrs=i,yrs=i+1,bndF=c(0.9,1.1))$hvt)
simHCR[["bound F"]]=bd

plot(simHCR)


## ----,fig.margin=TRUE----------------------------------------------------
bd=window(bd,end=30)
for (i in seq(29,49,1))
  bd=fwd(bd,catch=hcr(bd,refYrs=i,yrs=i+1,tac=T,bndTac=c(0.9,1.1))$tac)
simHCR[["bound TAC"]]=bd

plot(simHCR)


## ----hcr-ts,fig.width=6,fig.height=8-------------------------------------
p.=plot(simHCR[c(1,3)])+
  theme(legend.position="bottom")+
  scale_colour_manual(values=c("cyan4","grey10"),
                      labels=c("HCR","10% Constraint on Effort"))+
  guides(col=guide_legend(title=NULL))+
  scale_x_continuous(limits=c(27,50),breaks=c(30,40,50),labels=c(1,11,21))+
  scale_y_continuous(breaks=NULL)

p.$data=transform(p.$data,qname=factor(qname,levels=c("Harvest","Yield","Stock")))
p.+geom_hline(aes(yintercept=X1),data=
                cbind("qname"=c("Yield","Harvest","Stock"),
                      data.frame(refpts(bd))),col="red")


## ----hcrI,fig.width=5,fig.height=5---------------------------------------

kb=ldply(simHCR[c(1)],function(x,sd=.1)
  model.frame(mcf(FLQuants(stock  =rlnorm(100,log(stock(  x)%/%bmsy(x)),sd),
                           harvest=rlnorm(100,log(harvest(x)%/%fmsy(x)),sd)))))
kb=subset(kb,year%in% 29:50)

pt=ldply(simHCR[c(1)],function(x)
  model.frame(mcf(FLQuants(stock  =stock(  x)%/%bmsy(x),
                           harvest=harvest(x)%/%fmsy(x)))))
pt =subset(pt,year%in% 1:50)
pt.=ddply(pt,.(year,.id),with,data.frame(stock=median(stock),harvest=median(harvest)))

i=40
print(kobePhase(subset(kb,year%in%i))+
        geom_line(aes(stock,harvest),data=biodyn:::hcrPlot(bd),col="brown",size=1.5)+
        ggtitle(paste("Year",i-29)) +
        theme(legend.position="none",plot.background=element_rect(fill="transparent",colour=NA),
              plot.title = element_text(lineheight=.8, face="italic"))+
        xlim(0,2)+ylim(0,2)+
        geom_point(aes(stock,harvest,col=.id,fill=.id), shape=21, size = 3) + 
        geom_point(aes(stock,harvest,group=.id,fill=.id),col="black", data=subset(pt.,year==i),shape=21,size=5)+
        geom_path( aes(stock,harvest,col=.id,group=.id) ,data=subset(pt, year<=i),size=1)+
        scale_fill_manual(  values=c("cyan1","green","red","yellow")) + 
        scale_colour_manual(values=c("cyan4","green","red","yellow"))+
        coord_cartesian(xlim=c(0,2),ylim=c(0,2)))


## ----hcrII,fig.width=5,fig.height=5,cache=FALSE--------------------------
library(FLCore)
library(plyr)
library(kobe)
library(ggplot2)

kb=ldply(simHCR[1:2],function(x,sd=.1)
  model.frame(mcf(FLQuants(stock  =stock(  x)%/%bmsy(x),
                           harvest=harvest(x)%/%fmsy(x)))))
kb=subset(kb,year%in% 29:50)

pt=ldply(simHCR,function(x)
  model.frame(mcf(FLQuants(stock  =stock(  x)%/%bmsy(x),
                           harvest=harvest(x)%/%fmsy(x)))))
pt =subset(pt,year%in% 1:50)
pt.=ddply(pt,.(year,.id),with,data.frame(stock=median(stock),harvest=median(harvest)))

i=40
kobePhase(subset(kb,year%in%i))+
  # geom_line(aes(stock,harvest),data=hcrPlot(bd),col="brown",size=1.5)+
  ggtitle(paste("Year",i-29)) +
  theme(legend.position="none",plot.background=element_rect(fill="transparent",colour=NA),
        plot.title = element_text(lineheight=.8, face="italic"))+
  xlim(0,2)+ylim(0,2)+
  geom_point(aes(stock,harvest,col=.id,fill=.id,size=.id), shape=21) + 
  geom_point(aes(stock,harvest,fill=.id,group=.id),data=subset(pt.,year==i),shape=21,col="black",size=5)+
  geom_path( aes(stock,harvest,col=.id,group=.id) ,data=subset(pt, year<=i),size=1)+
  scale_fill_manual(  values=c("cyan1","cyan1","grey90","green","red","yellow")) + 
  scale_size_manual(  values=c(1,3)) + 
  scale_colour_manual(values=c("cyan4","cyan4","grey20","green","red","yellow"))+
  coord_cartesian(xlim=c(0,2),ylim=c(0,2))


## ----hcr-ts-mc,fig.width=8,fig.height=8----------------------------------
pe=rlnorm(100,FLQuant(0,dimnames=list(year=1:50)),0.5)
simHCR[[1]]=fwd(simHCR[[1]],harvest=harvest(simHCR[[1]])[,ac(1:50)],pe=pe)
simHCR[[3]]=fwd(simHCR[[3]],harvest=harvest(simHCR[[3]])[,ac(1:50)],pe=pe)
p.=plot(simHCR[c(1,3)])+
  theme(legend.position="bottom")+
  scale_colour_manual(values=c("cyan4","grey10"),
                      labels=c("HCR","10% Constraint on Effort"))+
  guides(col=guide_legend(title=NULL))+
  scale_x_continuous(limits=c(27,50),breaks=c(30,40,50),labels=c(1,11,21))+
  scale_y_continuous(breaks=NULL)

p.$data=transform(p.$data,qname=factor(qname,levels=c("Harvest","Yield","Stock")))
p.+geom_hline(aes(yintercept=X1),data=
                cbind("qname"=c("Yield","Harvest","Stock"),data.frame(refpts(bd))),col="red")


## ----,eval=FALSE---------------------------------------------------------
## biodyn:::mseBiodyn


## ----, fig.width = 10, fig.height = 2, fig.fullwidth = TRUE, fig.cap = "Full width figure"----
qplot(wt, mpg, data=mtcars, colour=factor(cyl))


## ----, fig.cap = "Another figure"----------------------------------------
qplot(factor(cyl), mpg, data = mtcars, geom = "boxplot")


## ----,eval=FALSE---------------------------------------------------------
## biodyn:::mseBiodyn


## ----, fig.width = 10, fig.height = 2, fig.fullwidth = TRUE, fig.cap = "Full width figure"----
qplot(wt, mpg, data=mtcars, colour=factor(cyl))


## ----, fig.cap = "Another figure"----------------------------------------
qplot(factor(cyl), mpg, data = mtcars, geom = "boxplot")


## ----, results='asis'----------------------------------------------------
library(xtable)
options(xtable.comment = FALSE)
options(xtable.booktabs = TRUE)
xtable(head(mtcars[,1:6]), caption = "First rows of mtcars")


