  #include "admodel.h"
  #include <string>
  #include <dnorm.cpp> //include functions from custom library
  using std::string;
  const double pi = 3.141592654;
  int mcmc_iteration = 0;
  int phz;    // phase
  double lb;  // lower bound
  double ub;  // upper bound
  ofstream mcmc_par("mcmc_par.csv");
  ofstream mcmc_bio("mcmc_bio.csv");
  ofstream priors("prr.txt");
  ofstream bounds("bnd.txt");
  ofstream    lls("lls.txt");
#include <admodel.h>
#include <contrib.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <pella.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  nc.allocate("nc");
  Cdata.allocate(1,2,1,nc,"Cdata");
  ni.allocate("ni");
  nIdx.allocate("nIdx");
  Idata.allocate(1,3,1,ni,"Idata");
  Cyear.allocate(1,nc);
  C.allocate(1,nc);
  Iyear.allocate(1,ni);
  Idx.allocate(1,ni);
  I.allocate(1,ni);
  X.allocate(1,ni);
  logI.allocate(1,ni);
 string run_name = string(adprogram_name);
 if(option_match(argc,argv,"-ind") > -1){
   run_name = argv[option_match(argc,argv,"-ind") + 1];
   run_name = run_name.substr(0, run_name.rfind("."));}
 change_datafile_name((adstring)run_name.c_str() + ".ctl");
  _r_plui.allocate(1,4,"_r_plui");
  _k_plui.allocate(1,4,"_k_plui");
  _p_plui.allocate(1,4,"_p_plui");
  _a_plui.allocate(1,4,"_a_plui");
  qPh.allocate(1,nIdx,"qPh");
  qLo.allocate(1,nIdx,"qLo");
  qHi.allocate(1,nIdx,"qHi");
  qPr.allocate(1,nIdx,"qPr");
  sPh.allocate(1,nIdx,"sPh");
  sLo.allocate(1,nIdx,"sLo");
  sHi.allocate(1,nIdx,"sHi");
  sPr.allocate(1,nIdx,"sPr");
 change_datafile_name((adstring)run_name.c_str() + ".prr");
  r_prior.allocate(1,4,"r_prior");
  k_prior.allocate(1,4,"k_prior");
  p_prior.allocate(1,4,"p_prior");
  a_prior.allocate(1,4,"a_prior");
  msy_prior.allocate(1,4,"msy_prior");
  bmsy_prior.allocate(1,4,"bmsy_prior");
  fmsy_prior.allocate(1,4,"fmsy_prior");
  q_prior.allocate(1,4,"q_prior");
  s_prior.allocate(1,4,"s_prior");
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
 phz = (int) _r_plui[1];
 lb  =       _r_plui[2];
 ub  =       _r_plui[3];
  _r.allocate(lb,ub,phz,"_r");
 phz = (int) _k_plui[1];
 lb  =       _k_plui[2];
 ub  =       _k_plui[3];
  _k.allocate(lb,ub,phz,"_k");
 phz = (int) _a_plui[1];
 lb  =       _a_plui[2];
 ub  =       _a_plui[3];
  _a.allocate(lb,ub,phz,"_a");
 phz = (int) _p_plui[1];
 lb  =       _p_plui[2];
 ub  =       _p_plui[3];
  _p.allocate(lb,ub,phz,"_p");
  _q.allocate(1,nIdx,qLo,qHi,qPh,"_q");
  _s.allocate(1,nIdx,sLo,sHi,sPh,"_s");
  r.allocate("r");
  k.allocate("k");
  a.allocate("a");
  p.allocate("p");
  q.allocate(1,nIdx,"q");
  s.allocate(1,nIdx,"s");
  cnow.allocate("cnow");
  bnow.allocate("bnow");
  fnow.allocate("fnow");
  msy.allocate("msy");
  bmsy.allocate("bmsy");
  fmsy.allocate("fmsy");
  cmsy.allocate("cmsy");
  bbmsy.allocate("bbmsy");
  ffmsy.allocate("ffmsy");
  bk.allocate("bk");
  fr.allocate("fr");
  bratio.allocate("bratio");
  fratio.allocate("fratio");
  slopeb.allocate("slopeb");
  slopef.allocate("slopef");
  ll.allocate(1,nIdx,"ll");
  SS.allocate("SS");
  #ifndef NO_AD_INITIALIZE
  SS.initialize();
  #endif
  s_.allocate("s_");
  #ifndef NO_AD_INITIALIZE
  s_.initialize();
  #endif
  nii.allocate(1,2,"nii");
  #ifndef NO_AD_INITIALIZE
    nii.initialize();
  #endif
  _bratio.allocate("_bratio");
  #ifndef NO_AD_INITIALIZE
  _bratio.initialize();
  #endif
  _fratio.allocate("_fratio");
  #ifndef NO_AD_INITIALIZE
  _fratio.initialize();
  #endif
  xy.allocate("xy");
  #ifndef NO_AD_INITIALIZE
  xy.initialize();
  #endif
  x.allocate("x");
  #ifndef NO_AD_INITIALIZE
  x.initialize();
  #endif
  y.allocate("y");
  #ifndef NO_AD_INITIALIZE
  y.initialize();
  #endif
  xx.allocate("xx");
  #ifndef NO_AD_INITIALIZE
  xx.initialize();
  #endif
  B.allocate(1,nc+1,"B");
  #ifndef NO_AD_INITIALIZE
    B.initialize();
  #endif
  Bfit.allocate(1,ni,"Bfit");
  #ifndef NO_AD_INITIALIZE
    Bfit.initialize();
  #endif
  Ifit.allocate(1,ni,"Ifit");
  #ifndef NO_AD_INITIALIZE
    Ifit.initialize();
  #endif
  RSS.allocate(1,nIdx,"RSS");
  #ifndef NO_AD_INITIALIZE
    RSS.initialize();
  #endif
  summary.allocate(1,nc,1,7,"summary");
  #ifndef NO_AD_INITIALIZE
    summary.initialize();
  #endif
  lpbnow.allocate("lpbnow");
  lpfnow.allocate("lpfnow");
  lpmsy.allocate("lpmsy");
  lpbmsy.allocate("lpbmsy");
  lpfmsy.allocate("lpfmsy");
  lpcmsy.allocate("lpcmsy");
  lpbbmsy.allocate("lpbbmsy");
  lpffmsy.allocate("lpffmsy");
  lpbratio.allocate("lpbratio");
  lpfratio.allocate("lpfratio");
  lpslopeb.allocate("lpslopeb");
  lpslopef.allocate("lpslopef");
  lpr.allocate("lpr");
  lpk.allocate("lpk");
  lpfr.allocate("lpfr");
  lpbk.allocate("lpbk");
  neglogL.allocate("neglogL");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

  admaster_slave_variable_interface(*this);
  halfnlog2pi = 0.5*ni*log(2*pi);
  nReg=5;
  stepN =50;
  stepSz=0.05;
  // Data
  Cyear = (ivector) row(Cdata,1);
  C     =           row(Cdata,2);
  
  Iyear = (ivector) row(Idata,1);
  I     =           row(Idata,2);
  Idx   = (ivector) row(Idata,3);
  logI  = log(I);
  X     = Iyear - Cyear[1] + 1;
  
  // Parameters
  _r    = _r_plui[4];
  _k    = _k_plui[4];
  _a    = _a_plui[4];
  _p    = _p_plui[4];
  for (int j=1; j<=nIdx; j++){
    // change from log
    //logq(j) = qPr[j];
    //logs(j) = sPr[j];
    _q(j) = qPr[j];
    _s(j) = sPr[j];
    }
  // likelihood profile 
  lpbnow  =bnow;
  lpfnow  =fnow;
  if (1!=2){
  lpmsy   =msy;
  lpbmsy  =bmsy;
  lpfmsy  =fmsy;
  lpcmsy  =cmsy;
  lpbbmsy =bbmsy;
  lpffmsy =ffmsy;
  lpbratio=bratio;
  lpfratio=fratio;
  lpslopeb=slopeb;
  lpslopef=slopef;
  lpr     =_r;
  lpk     =_k;
  lpfr    =fr;
  lpbk    =bk;
  }
  lpbnow.set_stepnumber(stepN);
  lpbnow.set_stepsize(stepSz);
  lpfnow.set_stepnumber(stepN);
  lpfnow.set_stepsize(stepSz);
  if (1!=2){
  lpmsy.set_stepnumber(stepN);
  lpmsy.set_stepsize(stepSz);
  lpfmsy.set_stepnumber(stepN);
  lpfmsy.set_stepsize(stepSz);
  lpbmsy.set_stepnumber(stepN);
  lpbmsy.set_stepsize(stepSz);
  lpcmsy.set_stepnumber(stepN);
  lpcmsy.set_stepsize(stepSz);
  lpbbmsy.set_stepnumber(stepN);
  lpbbmsy.set_stepsize(stepSz);
  lpffmsy.set_stepnumber(stepN);
  lpffmsy.set_stepsize(stepSz);
  lpbratio.set_stepnumber(stepN);
  lpbratio.set_stepsize(stepSz);
  lpfratio.set_stepnumber(stepN);
  lpfratio.set_stepsize(stepSz);
  lpslopeb.set_stepnumber(stepN);
  lpslopeb.set_stepsize(stepSz);
  lpslopef.set_stepnumber(stepN);
  lpslopef.set_stepsize(stepSz);
  lpr.set_stepnumber(stepN);
  lpr.set_stepsize(stepSz);
  lpk.set_stepnumber(stepN);
  lpk.set_stepsize(stepSz);
  lpfr.set_stepnumber(stepN);
  lpfr.set_stepsize(stepSz);
  lpbk.set_stepnumber(stepN);
  lpbk.set_stepsize(stepSz);
  }
 
}

void model_parameters::userfunction(void)
{
  neglogL =0.0;
  get_fit();
  get_neglogL();
  if(mceval_phase())
    write_mcmc();
  // likelihood profile 
  lpbnow  =bnow;
  lpfnow  =fnow;
  if (1!=2){
  lpmsy   =msy;
  lpbmsy  =bmsy;
  lpfmsy  =fmsy;
  lpcmsy  =cmsy;
  lpbbmsy =bbmsy;
  lpffmsy =ffmsy;
  lpbratio=bratio;
  lpfratio=fratio;
  lpslopeb=slopeb;
  lpslopef=slopef;
  lpr     =_r;
  lpk     =_k;
  lpfr    =fr;
  lpbk    =bk;
  }
}

void model_parameters::report()
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  write_bounds();
  write_priors();
  summary.initialize();
  get_summary();
  report<<setprecision(12)
        <<"# r"      <<endl<<r      <<endl
        <<"# k"      <<endl<<k      <<endl
        <<"# b0"     <<endl<<a      <<endl
        <<"# p"      <<endl<<p      <<endl
        <<"# q"      <<endl<<q      <<endl
        <<"# s"      <<endl<<s      <<endl
        <<"# RSS"    <<endl<<RSS    <<endl
        <<"# neglogL"<<endl<<neglogL<<endl<<endl;
  report<<setprecision(12)
        <<"# Model summary"<<endl
        <<" year stock catch index hat stockHat stock."<<endl
        <<summary<<endl;
 write_ll();
}

void model_parameters::get_fit(void)
{
  //r = mfexp(logr);
  r = _r;
  k = _k;
  a = _a;
  p = _p;
  for (int j=1; j<=nIdx; j++){
    //change from log
    //q[j]   = mfexp(logq[j]);
    //s[j]   = mfexp(logs[j]);
    q[j]   = _q[j];
    s[j]   = _s[j];
    }
  B[1] = k*a;
  for(int t=1; t<=nc; t++)
    B[t+1] = sfabs(B[t] + r/p*B[t]*(1-pow(B[t]/k,p)) - C[t]);
  for (int j=1; j<=ni; j++)
     //Ifit[j] = B(X[j])*q(Idx[j]);
     Ifit[j] = 0.5*(B(X[j])+B(X[j]+1))*q(Idx[j]);
  cnow =C[nc];
  fnow =C[nc]/B[nc];
  bnow =B[nc];
  msy  =r*k*pow(1/(1+p),1/p+1);
  bmsy =(k*pow((1/(1+p)),(1/p)));
  fmsy =msy/bmsy;
  cmsy =C[nc]	/msy;
  bbmsy=bnow/bmsy;
  ffmsy=fnow/fmsy;
  bk=bnow/k;
  fr=fnow/r;
   bratio=0.0;
   fratio=0.0;
  _bratio=0.0;
  _fratio=0.0;
  for (int i=nc; i>nc-3; i--){
    bratio+=B[i];
    fratio+=C[i]/B[i];
    _bratio+=B[i-3];
    _fratio+=C[i-3]/B[i-3];
    }
 bratio=bratio/_bratio;
 fratio=fratio/_fratio;
  xy=0.0; 
  x =0.0;
   y=0.0; 
  xx=0.0; 
  for (int i=nc; i>nc-nReg; i--){
    x +=i;
    xx+=i*i;  
    y +=B[i];  
    xy+=i*B[i];  
    }
  slopeb = (nReg*xy - x*y)/(nReg*xx - x*2.0);
  xy=0.0; 
  x =0.0;
   y=0.0; 
  xx=0.0; 
  for (int i=nc; i>nc-nReg; i--){
    x +=i;
    xx+=i*i;  
    y +=C[i]/B[i];  
    xy+=i*C[i]/B[i];  
    }
  slopef = (nReg*xy - x*y)/(nReg*xx - x*2.0);
}

void model_parameters::get_neglogL(void)
{
  neglogL = halfnlog2pi;
  for (int j=1; j<=ni; j++){	
     //change from log
     //s_       = mfexp(logs[Idx[j]]);
     s_       = _s[Idx[j]];
     neglogL += log(s_)   
             +  pow(log(I[j])-log(Ifit[j]),2.0)/(2*s_*s_);
     }
  //neglogL = halfnlog2pi + ni*log(s(1)) + RSS[1]/(2*s(1)*s(1));
  // weighted likelihood priors
  dvariable _msy  =r*k*pow(1/(1+p),1/p+1);
  dvariable _bmsy =(k*pow((1/(1+p)),(1/p)));
  dvariable _fmsy =msy/bmsy;
  if (r_prior[1]>0)    neglogL += r_prior[1]*dnorm(r, r_prior[2], r_prior[3]); // /dnorm(r_prior[2], r_prior[2], r_prior[3]);
  if (k_prior[1]>0)    neglogL += k_prior[1]*dnorm(k, k_prior[2], k_prior[3]); // /dnorm(k_prior[2], k_prior[2], k_prior[3]);
  if (p_prior[1]>0)    neglogL += p_prior[1]*dnorm(p, p_prior[2], p_prior[ 3]); // /dnorm(p_prior[2], p_prior[2], p_prior[3]);
  if (a_prior[1]>0)    neglogL += a_prior[1]*dnorm(a, a_prior[2], a_prior[3]); // /dnorm(a_prior[2], a_prior[2], a_prior[3]);
  if ( msy_prior[1]>0) neglogL += msy_prior[1]*dnorm(_msy,  msy_prior[2],  msy_prior[3]); // /dnorm( msy_prior[2],  msy_prior[2],  msy_prior[3]);
  if (bmsy_prior[1]>0) neglogL += bmsy_prior[1]*dnorm(_bmsy, bmsy_prior[2], bmsy_prior[3]); // /dnorm(bmsy_prior[2], bmsy_prior[2], bmsy_prior[3]);
  if (fmsy_prior[1]>0) neglogL += fmsy_prior[1]*dnorm(_fmsy, fmsy_prior[2], fmsy_prior[3]); // /dnorm(fmsy_prior[2], fmsy_prior[2], fmsy_prior[3]);
 //for (i=1; i<=nIdx; i++){
  //  if (q_prior[i]>0) neglogL += q_prior[1]*dnorm(q, q_prior[2], q_prior[3]); // /dnorm(q_prior[2], q_prior[2], q_prior[3]);
  //  if (s_prior[i]>0) neglogL += s_prior[1]*dnorm(s, s_prior[2], s_prior[3]); // /dnorm(s_prior[2], s_prior[2], s_prior[3]);}
}

void model_parameters::get_summary(void)
{
  summary.colfill(1,(dvector)Cyear);
  summary.colfill(2,B);
  summary.colfill(3,C);
  for(int i=1; i<=ni; i++)  // allow missing years in abundance index
    {
    summary(X[i],4) = I[i];
    summary(X[i],5) = Ifit[i];
    summary(X[i],6) = Ifit[i]/q(Idx(i));
    summary(X[i],7) = (B[X[i]]+B[X[i]+1])/2.0;
    }
}

void model_parameters::write_mcmc(void)
{
  mcmc_iteration++;
  // Parameters
  if(mcmc_iteration == 1){
    mcmc_par<<"neglogL,r,k,b0,p,q,s"<<endl;
    mcmc_par<<neglogL<<","
          <<r      <<","
          <<k      <<","
          <<a      <<","
          <<p      <<","
          <<q      <<","
          <<s  <<endl;}
  // Biomass
  if(mcmc_iteration == 1){
    mcmc_bio<<Cyear[1];
    for(int t=2; t<=nc; t++)
      mcmc_bio<<","<<Cyear[t];
     mcmc_bio<<endl;
     }
  mcmc_bio<<B[1];
  for(int t=2; t<=nc; t++)
    mcmc_bio<<","<<B[t];
  mcmc_bio<<endl;
}

void model_parameters::write_priors(void)
{
  priors<<r_prior<<endl;
}

void model_parameters::write_bounds(void)
{
  bounds<<q	
        <<s<<endl;	
}

void model_parameters::write_ll(void)
{
   for (int j=1; j<=nIdx; j++) 
       ll[j]=0.0; 
    for (int j=1; j<=ni; j++){
       s_    = _s[Idx[j]];
       ll[Idx[j]] += log(s_) +pow(log(I[j])-log(Ifit[j]),2.0)/(2*s_*s_);
       }
  lls<<"# ll"<<endl<<ll<<endl;
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

void model_parameters::set_runtime(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  arrmblsize = 40000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000);
  gradient_structure::set_MAX_NVAR_OFFSET(10000);
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
