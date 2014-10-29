
## ----initKnitr, echo=FALSE, results="asis"-------------------------------
library(knitr)

#output: rmarkdown::tufte_handout

# output:
#   rmdformats::html_clean:
#     fig_width: 6
#     fig_height: 6
#     highlight: pygments

## Global options
options(max.print="75")
opts_chunk$set(fig.path="out/",
               echo =TRUE,
               eval=TRUE,
               cache=TRUE,
               cache.path="cache/",
               prompt=FALSE,
               tidy=FALSE,
               comment=NA,
               message=FALSE,
               warning=FALSE,
               fig.margin=TRUE, 
               fig.height=6, 
               fig.width=4)

opts_knit$set(width=75)


## ----,fig.caption="Simulated stock"--------------------------------------
library(plyr)
library(biodyn)

bd=sim()
bd=window(bd,end=49)


## ----, fig,caption="Simulated CPUE series",fig.height=2------------------
cpue=(stock(bd)[,-dims(bd)$year]+stock(bd)[,-1])/2
cpue=rlnorm(1,log(cpue),.2)

ggplot(cpue)+geom_point(aes(year,data))


## ----, fig.cap="Fitted compared to actual time series"-------------------
#set parameters
setParams(bd) =cpue
setControl(bd)=params(bd)
control(bd)[3:4,"phase"]=-1


## ----, fig.cap="Fitted compared to actual time series"-------------------
#fit
bdHat=fit(bd,cpue)

plot(biodyns("Estimate"=bdHat,"Actual"=bd))+
  theme(legend.position="bottom")


## ----, fig.cap="Jack knife error bars"-----------------------------------
bdJK=fit(bdHat,jackknife(cpue))


## ----, fig.cap="Jack knife error bars"-----------------------------------
plotJack(bdHat,bdJK)


## ------------------------------------------------------------------------
#FLQuant
library(ggplotFL)

true=stock(bd   )[,49]
hat =stock(bdHat)[,49]
jack=stock(bdJK )[,49]


n   =dims(jack)$iter
mn  =apply(jack, 1:5, mean)   
rsdl=sweep(jack,1:5,mn,"-")
ss  =apply(rsdl^2,1:5,sum)
            
bias         =(n-1)*(hat-mn)
biasCorrected=n*hat-(n-1)*mn

se  =sqrt(((n-1)/n)*ss)


## ------------------------------------------------------------------------
# FLPar
true=refpts(bd   )
hat =refpts(bdHat)
jack=refpts(bdJK )

n   =dims(jack)$iter
mn  =apply(jack,  seq(length(dim(jack))-1), mean)   
rsdl=sweep(jack,  seq(length(dim(jack))-1), mn,"-")
ss  =apply(rsdl^2,seq(length(dim(jack))-1),sum)
            
bias         =(n-1)*(hat-mn)
biasCorrected=n*hat-(n-1)*mn

se  =sqrt(((n-1)/n)*ss)


