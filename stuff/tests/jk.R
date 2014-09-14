library(biodyn)
library(reshape)
library(plyr)

source('~/Desktop/flr/git/biodyn/R/biodyn-jackRand.R')
source('~/Desktop/flr/git/biodyn/R/biodyn-jackSummary.R')

set.seed(987)
bd  =window(simBiodyn(),end=40)
cpue=rlnorm(1,log((stock(bd)[,-40]+stock(bd)[,-1])/2),0.2)

setParams( bd)=cpue
setControl(bd)=params(bd)

bd@control[3:4,"phase"]=-1

bdHat=fit(bd,cpue)

plot(biodyns("True"=bd,"Hat"=bdHat))

bdJk=fit(bdHat,jackknife(cpue))

plot(biodyns("True"=bd,"Hat"=bdHat,"Jack"=bdJk))

jackSummary(params(bdHat),params(bdJk))$cov

bdMC=randJack(100,bdHat,bdJk)

plot(biodyns("MC"=bdMC,"Jack"=bdJk))

