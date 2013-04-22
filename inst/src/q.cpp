
calcSigma=function(stock,idx=rep(0,length(obs)),error="log",na.rm=T){
      if (error=="log") stock=log(stock)
      if (error=="log") idx  =log(idx)

      ####  Biomass mid year
      mnBio=function(x) (x[-length(x)]+x[-1])/2
          
      ####  Biomass mid year
      stock=mnBio(stock)
      idx  =idx[seq(length(stock))]
      hat  =predict(lm(idx~stock))
              
      SS =sum((idx[as.numeric(names(hat))]-hat)^2,na.rm=na.rm)
          
      return((SS/length(hat))^.5)}
       

double calcQ(double *stock, double *index, int n, error=1) { 

  double q;

  \\ Biomass mid year
  for (i=0; i<n; i++)
    stock[i]=(stock[i]+stock[i+1])/2.0;
  
  \\normal 
  if (error==1){
      double si =0.0, ii=0.0;
      for (i=0; i<n; i++){
        si+=stock[i]*index[i];
        ii+=index[i]*index[i];
        }
     q    =si/ii;}
  \\lognormal
  else if (error==2){
    double si=0.0;
    for (i=0; i<n; i++)
        si+=log(stock)-log(index);    
     
    q     =exp(si)/n;}
  \\cv
  else if (error==3){
    res   =sum(hat/obs)
    sigma2=calcSigma(res,na.rm=na.rm)
    q     =(-res+(res^2+4*n*sigma2*sum((hat/obs)^2)))/(2*n*sigma2)
  }
  
  return q;}




