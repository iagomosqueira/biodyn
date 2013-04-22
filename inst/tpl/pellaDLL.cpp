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
#include <admodel.h>

  extern "C"  {
    void ad_boundf(int i);
  }
#include <pella.htp>

model_data::model_data(int argc,char * argv[],dll_args& ad_dll) : ad_comm(argc,argv)
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
  logr_plui.allocate(1,4,"logr_plui");
  logk_plui.allocate(1,4,"logk_plui");
  logp_plui.allocate(1,4,"logp_plui");
  loga_plui.allocate(1,4,"loga_plui");
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
  q_prior.allocate(1,4,"q_prior");
  s_prior.allocate(1,4,"s_prior");
}

model_parameters::model_parameters(int sz,int argc,char * argv[], dll_args& ad_dll) : 
 model_data(argc,argv,ad_dll) , function_minimizer(sz)
{
  initializationfunction();
 phz = (int) logr_plui[1];
 lb  =       logr_plui[2];
 ub  =       logr_plui[3];
  logr.allocate(lb,ub,phz,"logr");
 phz = (int) logk_plui[1];
 lb  =       logk_plui[2];
 ub  =       logk_plui[3];
  logk.allocate(lb,ub,phz,"logk");
 phz = (int) loga_plui[1];
 lb  =       loga_plui[2];
 ub  =       loga_plui[3];
  loga.allocate(lb,ub,phz,"loga");
 phz = (int) logp_plui[1];
 lb  =       logp_plui[2];
 ub  =       logp_plui[3];
  logp.allocate(lb,ub,phz,"logp");
  logq.allocate(1,nIdx,qLo,qHi,qPh,"logq");
  logs.allocate(1,nIdx,sLo,sHi,sPh,"logs");
  r.allocate("r");
  k.allocate("k");
  a.allocate("a");
  p.allocate("p");
  q.allocate(1,nIdx,"q");
  s.allocate(1,nIdx,"s");
  Ynow.allocate("Ynow");
  Bnow.allocate("Bnow");
  Hnow.allocate("Hnow");
  MSY.allocate("MSY");
  BMSY.allocate("BMSY");
  HMSY.allocate("HMSY");
  YMSY.allocate("YMSY");
  BBMSY.allocate("BBMSY");
  HHMSY.allocate("HHMSY");
  Bk.allocate("Bk");
  Hr.allocate("Hr");
  BRatio.allocate("BRatio");
  HRatio.allocate("HRatio");
  slopeB.allocate("slopeB");
  slopeH.allocate("slopeH");
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
  _BRatio.allocate("_BRatio");
  #ifndef NO_AD_INITIALIZE
  _BRatio.initialize();
  #endif
  _HRatio.allocate("_HRatio");
  #ifndef NO_AD_INITIALIZE
  _HRatio.initialize();
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
  lpr.allocate("lpr");
  neglogL.allocate("neglogL");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
}

void model_parameters::preliminary_calculations(void)
{

  admaster_slave_variable_interface(*this);
  halfnlog2pi = 0.5*ni*log(2*pi);
  nReg=5;
  // Data
  Cyear = (ivector) row(Cdata,1);
  C     =           row(Cdata,2);
  
  Iyear = (ivector) row(Idata,1);
  I     =           row(Idata,2);
  Idx   = (ivector) row(Idata,3);
  logI  = log(I);
  X     = Iyear - Cyear[1] + 1;
  
  // Parameters
  logr    = logr_plui[4];
  logk    = logk_plui[4];
  loga    = loga_plui[4];
  logp    = logp_plui[4];
  for (int j=1; j<=nIdx; j++){
    logq(j) = qPr[j];
    logs(j) = sPr[j];}
  // likelihood profile 
  lpr   =mfexp(logr);
  
}

void model_parameters::userfunction(void)
{
  neglogL =0.0;
  get_fit();
  get_neglogL();
  if(mceval_phase())
    write_mcmc();
  // likelihood profile 
  lpr=mfexp(logr);
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
}

void model_parameters::get_fit(void)
{
  r = mfexp(logr);
  k = mfexp(logk);
  a = mfexp(loga);
  p = mfexp(logp);
  for (int j=1; j<=nIdx; j++){
    q[j]   = mfexp(logq[j]);
    s[j]   = mfexp(logs[j]);
    }
  B[1] = k*a;
  for(int t=1; t<=nc; t++)
    B[t+1] = sfabs(B[t] + r/p*B[t]*(1-pow(B[t]/k,p)) - C[t]);
  for (int j=1; j<=ni; j++)
     //Ifit[j] = B(X[j])*q(Idx[j]);
     Ifit[j] = 0.5*(B(X[j])+B(X[j]+1))*q(Idx[j]);
  Ynow =C[nc];
  Hnow =C[nc]/B[nc];
  Bnow =B[nc];
  MSY  =r*k*pow(1/(1+p),1/p+1);
  BMSY =(k*pow((1/(1+p)),(1/p)));
  HMSY =MSY/BMSY;
  YMSY =C[nc]	/MSY;
  BBMSY=Bnow/BMSY;
  HHMSY=Hnow/HMSY;
  Bk=Bnow/k;
  Hr=Hnow/r;
   BRatio=0.0;
   HRatio=0.0;
  _BRatio=0.0;
  _HRatio=0.0;
  for (int i=nc; i>nc-3; i--){
    BRatio+=B[i];
    HRatio+=C[i]/B[i];
    _BRatio+=B[i-3];
    _HRatio+=C[i-3]/B[i-3];
    }
 BRatio=BRatio/_BRatio;
 HRatio=HRatio/_HRatio;
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
  slopeB = (nReg*xy - x*y)/(nReg*xx - x*2.0);
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
  slopeH = (nReg*xy - x*y)/(nReg*xx - x*2.0);
}

void model_parameters::get_neglogL(void)
{
  neglogL = halfnlog2pi;
  for (int j=1; j<=ni; j++){	
     s_       = mfexp(logs[Idx[j]]);
     neglogL += log(s_)   
             +  pow(log(I[j])-log(Ifit[j]),2.0)/(2*s_*s_);
     }
  //neglogL = halfnlog2pi + ni*log(s(1)) + RSS[1]/(2*s(1)*s(1));
  // weighted likelihood priors
  if (r_prior[1]>1) neglogL -= dnorm(r, r_prior[2], r_prior[3])/dnorm(r_prior[2], r_prior[2], r_prior[3]);
  if (k_prior[1]>1) neglogL -= dnorm(k, k_prior[2], k_prior[3])/dnorm(k_prior[2], k_prior[2], k_prior[3]);
  if (p_prior[1]>1) neglogL -= dnorm(p, p_prior[2], p_prior[3])/dnorm(p_prior[2], p_prior[2], p_prior[3]);
  if (a_prior[1]>1) neglogL -= dnorm(a, a_prior[2], a_prior[3])/dnorm(a_prior[2], a_prior[2], a_prior[3]);
  //for (i=1; i<=nIdx; i++){
  //  neglogL += q_prior[i]*dnorm(q[i],    q_prior[i,2],     q_prior[i,3]);
  //  neglogL += s_prior[i]*dnorm(s[i],    s_prior[i,2],     s_prior[i,3]);}
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
  bounds<<logq
        <<logs
        <<q	
        <<s<<endl;	
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
extern "C" {

#if !defined(__MSVC32__)
#  define __declspec(x)
#endif

#if !defined(__BORLANDC__)
#  define _export
#else
#  define _export __stdcall
#endif

__declspec(dllexport) void _export pella(char ** dll_options)
{
  arrmblsize = 40000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000);
  gradient_structure::set_MAX_NVAR_OFFSET(10000);
  int argc=1;
  try {
    char **argv=parse_dll_options("pella",argc,*dll_options);
    do_dll_housekeeping(argc,argv);
    dll_args ad_dll;
    gradient_structure::set_NO_DERIVATIVES();
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv,ad_dll);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
    ad_make_code_reentrant();
  }
  catch (spdll_exception spe){ 
    if (ad_printf && spe.e) (*ad_printf) ("abnormal exit from newtest\n");
  }
}
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
