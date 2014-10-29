
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
               cache=FALSE,
               cache.path="cache/",
               prompt=FALSE,
               tidy=FALSE,
               comment=NA,
               message=FALSE,
               warning=FALSE)

opts_knit$set(width=75)


## ----,eval=FALSE---------------------------------------------------------
## library(knitr)
## purl(    "/home/laurie/Desktop/flr/git/biodyn/stuff/vignettes/biodyn.Rmd",
##      out="/home/laurie/Desktop/flr/git/biodyn/vignettes/biodyn.R")
## 


## ------------------------------------------------------------------------
library(biodyn)
library(ggplot2)
library(plyr)


## ------------------------------------------------------------------------
bd =biodyn()


## ------------------------------------------------------------------------
bd=biodyn(catch=FLQuant(100,dimnames=list(year=1990:2010)))


## ----,eval=FALSE---------------------------------------------------------
## data(ple4)
## bd =as(ple4,"biodyn")


## ----,eval=FALSE---------------------------------------------------------
## library(aspic)
## asp=aspic("http://http://rscloud.iccat.int/kobe/swon/2013/aspic/run2/")
## bd =as(asp,"biodyn")


## ------------------------------------------------------------------------
bd=sim()


## ----, fig.margin=TRUE, fig.width=4, fig.height=8, fig.cap="Simulated time series"----
bd=sim()
bd=window(bd,end=49)
plot(bd)


## ----, fig.margin=TRUE, eval=FALSE, fig.height=6, fig.width=4, echo=FALSE, fig.cap="Simulated CPUE series"----
## biodyn:::plotEql(bd)


## ----, fig.margin=TRUE, fig.cap="Simulated CPUE series"------------------
library(reshape)
x=sim()
plotPrd(x)+
  geom_path( aes(stock,catch),
             model.frame(FLQuants(x,"stock","catch")))+
  geom_point(aes(stock,catch),
             model.frame(FLQuants(x,"stock","catch")))


## ----, fig.margin=TRUE, eval=FALSE, fig.cap="Simulated CPUE series"------
## plotMSE()


## ----, fig.margin=TRUE, fig.height=6, fig.width=4,fig.cap="Monte Carlo simuation time series with confidence intervals and a worm (single simulation)."----
harvest=rlnorm(200,log(harvest(bd)[,-1]),.2)

bd=fwd(bd,harvest=harvest)

plot(bd,worm=3)+
  theme(legend.position="bottom")


## ----, echo=TRUE, fig.margin=TRUE, fig.height=4, fig.cap="Simulated stock"----
bd=sim()


## ----, fig.margin=TRUE, fig.cap="Simulated CPUE series"------------------
cpue=(stock(bd)[,-dims(bd)$year]+
      stock(bd)[,-1])/2
cpue=rlnorm(1,log(cpue),.2)

ggplot(as.data.frame(cpue))+
  geom_point(aes(year,data))+
  geom_line(aes(year,data),data=as.data.frame(stock(bd)),col="salmon")


## ----,eval=FALSE---------------------------------------------------------
## bd=biodyn(catch=catch(bd),msy=mean(catch(bd)))


## ----,fig.margin=TRUE,fig.width=4,fig.height=6---------------------------
setParams( bd)=cpue
params(bd)


## ----,fig.margin=TRUE,fig.width=4,fig.height=6---------------------------
setControl(bd)=params(bd)
  
bd@control


## ----, eval=TRUE, fig.margin=TRUE, fig.height=6,fig.cap="A comparison of the true and fitted time series"----
bd@control[3:4,"phase"]=-1
bdHat=fit(bd,cpue)

# plot(biodyns("True"=bd,"Hat"=bdHat))+
#   theme(legend.position="bottom")


## ----,fig.margin=TRUE,fig.width=4,fig.height=6---------------------------
params(bdHat)
params(bdHat)/params(bd)


## ----,echo=TRUE----------------------------------------------------------
head(bdHat@diags)


## -----rsdl5, fig.width=4,fig.height=4,fig.margin=TRUE, fig.cap="Quantile-quantile plot to compare residual distribution with the normal distribution."----
rsdl=bdHat@diags
ggplot(rsdl)                                           +
  geom_point( aes(qqx,qqy))                            +
  stat_smooth(aes(qqx,qqHat),method="lm",se=T,fill="blue", alpha=0.1)         +
  theme_ms(14,legend.position="bottom")               


## -----rsdl2, fig.margin=TRUE, fig.height=4, figwidth=4, fig.cap="Observed CPUE verses fitted, blue line is a linear resgression fitted to points, black the y=x line."----
ggplot(with(rsdl, data.frame(obs=stdz(index),hat=stdz(hat)))) +
    geom_abline(aes(0,1))                                     +
    geom_point( aes(obs,hat))                                 +
    stat_smooth(aes(obs,hat),method="lm", se=F)               +
    theme_ms(14,legend.position="bottom")                     +
    xlab("Fitted") + ylab("Observed")


## -----rsdl3,fig.height=3, fig.margin=TRUE, fig.cap="Residuals by year, with lowess smoother"----
dat=transform(subset(rsdl,!is.na(residual), 
                     residual=stdz(residual,na.rm=T)))

ggplot(aes(year,residual),data=dat)  +
  geom_hline(aes(yintercept=0))      +
  geom_point()                       +
  stat_smooth(method="loess",se=F)   +
  theme_ms(14,legend.position="bottom") 


## -----rsdl6,fig.height=3, fig.margin=TRUE, fig.cap="Plot of residuals against fitted value, to check variance relationship."----
ggplot(aes(hat, residual),
       data=subset(rsdl,!is.na(hat) & !is.na(residual)))   +
  geom_hline(aes(yintercept=0))         +
  geom_point()                          +
  stat_smooth(method="loess",se=F)      +
  theme_ms(14,legend.position="bottom")


## -----rsdl4, fig.width=4,fig.width=4,fig.margin=TRUE, fig.cap="Plot of autocorrelation, i.e. $residual_{t+1}$ verses $residual_{t}$."----
ggplot(rsdl)                                              +
  geom_point( aes(residual,residualLag))                  +
  stat_smooth(aes(residual,residualLag),method="lm",se=F) +
  geom_hline(aes(yintercept=0))     +
  xlab(expression(Residual[t]))     + 
  ylab(expression(Residual[t+1]))   +
  theme_ms(14,legend.position="bottom")  


## ----1d, fig.margin=TRUE, fig.cap="Likelihood profile for r"-------------
load("/home/laurie/Desktop/bdHat.RData")

bdHat=fit(bdHat,cpue)
setControl(bdHat)=params(bdHat)
res=profile(bdHat,which='r',fixed=c('b0','p'),
            cpue,range=seq(0.95,1.03,.002))
ggplot(subset(res,ll.u1<0))+geom_line(aes(r,ll.u1))


## ----2d, eval=FALSE, fig.margin=TRUE, fig.cap="Likelihood profile for r"----
## res=profile(bdHat,which=c('r','k'),fixed=c('b0','p'),
##             cpue,range=seq(0.97,1.03,.02))
## ggplot(res, aes(r, k, z=ll.u1))+
##   stat_contour(aes(colour = ..level..), size = 1)


## ----like, fig.margin=TRUE, fig.height=6, fig.width=4, fig.cap="Likelihood profile by data conmponent, i.e. CPUE series"----

bd=sim()

Us  =FLQuants("Unbiased"     =
                rlnorm(1,log((stock(bd)[,-dims(bd)$year]+
                              stock(bd)[,-1])/2),0.2),
              "Increase in q"=
                rlnorm(1,log((stock(bd)[,-dims(bd)$year]+
                              stock(bd)[,-1])/2),0.2))

bds=bd

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

ggplot(subset(melt(prfl[,c("r","ll.u1","ll.u2","1")],id="r"),
              value<10))+
  geom_line(aes(r,value,group=variable,col=variable))+
  facet_wrap(~variable,scale="free",ncol=1)          +
  theme(legend.position="bottom")


## ----,echo=FALSE,eval=FALSE----------------------------------------------
## bd2=fit(bdHat,cpue,cmdOps="-lprof")
## save(bd2,bdHat,cpue,file="/home/laurie/Desktop/bd2.RData")
## prf=subset(bd@profile, param %in% c("bbmsy","ffmsy"))
## prf=data.frame(What="Profile",t(daply(prf, .(param), with, sample(value,500,prob=p,replace=T))))
## names(prf)[2:3]=c("Stock","Harvest")


## ----, fig.margin=TRUE---------------------------------------------------
sims=biodyns("Best Fit"=bd)


## ----,fig.height=6, fig.margin=TRUE--------------------------------------
save(bdHat,cpue,file="/home/laurie/Desktop/bdHat.RData")
sims[["Vcov"]]=mvn(bdHat,500,nms=c("r","k"),fwd=TRUE)


## ----boot, fig.height=4, fig.margin=TRUE, fig.cap="Bootstrapped CPUE series"----
cpueSim=bdHat@diags[,c("year","hat")]
names(cpueSim)[2]="data"
cpueSim=as.FLQuant(cpueSim)

#cv(diags["residuals"])

cpueSim=rlnorm(100,log(cpueSim),0.25)

cpueSim[cpueSim==0]=NA

plot(cpueSim)

sims[["CPUE"]]=fit(propagate(bdHat,100),cpueSim)


## ----jackknife,fig.height=4,fig.margin=TRUE, fig.cap="Plot predicted stock trend by index"----
bdJK=fit(bdHat,jackknife(cpue))
save(bdJK,file="/home/laurie/Desktop/bdJK.RData")


## ----jackknife2,fig.height=4,fig.margin=TRUE, fig.cap="Predicted stock trend by index"----
source('~/Desktop/flr/git/biodyn/R/biodyn-jackRand.R')
sims[["Jack Knife"]]=randJack(500,bdJK)


## ------------------------------------------------------------------------
sims[["MCMC"]]=fit(bdHat,cpue,cmdOps=c("-mcmc 1000000, -mcsave 5000"))


## ----,fig.height=4,fig.margin=TRUE,--------------------------------------
acf(c(params(sims[["MCMC"]])["r"]))


## ----,eval=FALSE---------------------------------------------------------
## bdHat@mng


## ----,eval=FALSE---------------------------------------------------------
## bdHat@mngVcov


## ----,eval=FALSE,fig.margin=TRUE,fig.width=4,fig.height=6----------------
## currentState   =bdHat@mng[c("bbmsy","ffmsy"),"hat",drop=T]
## currentStateVar=bdHat@mngVcov[c("bbmsy","ffmsy"),
##                            c("bbmsy","ffmsy"),drop=T]
## 
## refs=mvrnorm(100,currentState,currentStateVar)
## 
## ggplot(data=as.data.frame(refs))+
##    geom_density(aes(x=bbmsy,y = ..count..))


## -----7,fig.fullwidth=TRUE, fig.width=4, fig.height=6,fig.cap="Densities of Stock/BMSY from different methods for estimating uncertainty."----
source('~/Desktop/flr/git/biodyn/R/biodyn-kobe.R')
df=kobe(sims[-1])
ggplot(subset(df,year==49))+ 
  geom_density(aes(x=stock, y=..count..), position = "stack",fill="red")+
  facet_wrap(~.id,scale="free_y",ncol=1)+
  geom_vline(aes(xintercept=c(stock(sims[[1]])[,"49"]%/%bmsy(sims[[1]]))))
  #scale_x_continuous(limits=c(0.5,1.5))


## -----8,fig.fullwidth=TRUE, fig.width=6, fig.height=8,fig.cap="Densities of Harvest/FMSY from different methods for estimating uncertainty."----

ggplot(subset(df,year==49))+ 
  geom_density(aes(x=harvest, y=..count..), position = "stack",fill="red")+
  facet_wrap(~.id,scale="free_y",ncol=1)+
  scale_x_continuous(limits=c(0,2.0))+
  geom_vline(aes(xintercept=c(harvest(sims[[1]])[,"48"]%/%fmsy(sims[[1]]))))


## -----9,fig.margin=TRUE,fig.width=4,fig.height=4,fig.caption="Kobe Phase Plots"----
kobePhase()+
  geom_point(aes(stock,harvest),data=subset(df,year==49))+
  facet_wrap(~.id,ncol=2)


## ----, fig.margin=TRUE,fig.width=4, fig.height=6,fig.cap="Projection"----
harvest=rlnorm(100,log(harvest(bdHat))[,-dims(bdHat)$year],.1)
bdHat =fwd(bdHat,harvest=harvest)

plot(bdHat,worm=c(2,8))+    
  theme(legend.position="bottom")


## ------------------------------------------------------------------------
bd   =sim()

bd=window(bd,end=29)
for (i in seq(29,49,1))
  bd=fwd(bd,harvest=hcr(bd,refYrs=i,yrs=i+1)$hvt)
simHCR=biodyns("1"=bd)


## ------------------------------------------------------------------------
bd=window(bd,end=29)
for (i in seq(29,49,3))
  bd=fwd(bd,harvest=hcr(bd,refYrs=i,yrs=i+1:3)$hvt)
simHCR[["3"]]=bd


## ------------------------------------------------------------------------
bd=window(bd,end=29)
for (i in seq(29,49,1))
  bd=fwd(bd,harvest=hcr(bd,refYrs=i,yrs=i+1,bndF=c(0.9,1.1))$hvt)
simHCR[["bound F"]]=bd


## ------------------------------------------------------------------------
bd=window(bd,end=30)
for (i in seq(29,49,1))
  bd=fwd(bd,catch=hcr(bd,refYrs=i,yrs=i+1,tac=T,bndTac=c(0.9,1.1))$tac)
simHCR[["bound TAC"]]=bd


## ----,fig.fullwidth=TRUE,fig.width=6,fig.height=6,fig.cap="Plots of projections"----
plot(simHCR)+
  theme(legend.position="bottom")


## -----hcrII,fig.margin=TRUE,fig.width=6,fig.height=6---------------------
pe=rlnorm(500,FLQuant(0,dimnames=list(year=1:50)),0.5)

bd=window(sim(),end=30)
bd.=bd
bd@stock =propagate(bd@stock, 500)
bd=fwd(bd,harvest=harvest(bd)[,2:30],pe=pe)

for (i in seq(30,48,1))
  bd=fwd(bd,
         catch=hcr(bd,refYrs=i,yrs=i+1,tac=T,bndTac=c(0.9,1.1))$tac,
         pe   =pe)

plot(bd)


## -----hcrKobe,fig.margin=TRUE,fig.width=6,fig.height=6-------------------
trks=biodyn::kobe(bd,what="trks")
trks=mdply(data.frame(Year=seq(33,49,3)), 
           function(Year) subset(trks,year<=Year))

pts =transform(biodyn::kobe(bd,what="pts",year=seq(33,49,3)),
                 Year=year)[,c("stock","harvest","Year")]

kobePhase()+    
    geom_line(aes(stock,harvest),data=hcrPlot(bd.),
              col="brown",size=1.5)                             +    
    geom_path( aes(stock,harvest),data=subset(trks,pctl=="50%"))+
    geom_point(aes(stock,harvest),data=subset(pts,Year>=33))    +
    facet_wrap(~Year)


## ----,eval=FALSE---------------------------------------------------------
## biodyn:::mseBiodyn


## ----, results='asis'----------------------------------------------------
library(xtable)
options(xtable.comment = FALSE)
options(xtable.booktabs = TRUE)
xtable(head(mtcars[,1:6]), caption = "First rows of mtcars")


