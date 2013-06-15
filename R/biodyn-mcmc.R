par2mcmc=function (x)
  mcmc(t(params(x)[drop=T]), start = 1, thin = attributes(params(x))$mcsave)
  