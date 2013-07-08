par2mcmc=function (x)
  mcmc(t(x[drop=T]), start = 1, thin = attributes(x)$mcsave)
  