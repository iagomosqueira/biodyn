//=========================================================================================================================
// File:        pella.tpl
// Model:       Pella-Tomlinson model, with Binit=k*a
// Parameters:  r, k, a, p, q, s
// Fitted data: Abundance index
// Likelihood:  Log-transformed normal
// References:  Polacheck et al. (1993)
// Notes:       q and s are free parameters, to allow full uncertainty in MCMC analysis
// History:      9 Mar 2010 Arni Magnusson created, to benchmark against R optimizers
//               7 Oct 2010 Arni Magnusson improved string handling and comments
//=========================================================================================================================
// Implementation notes
//   Abundance index may not exist for all years
//   Vectors that include all years: B, C
//   Vectors that include abundance index years: I, Ifit, X
//   X links long and short vectors
//=========================================================================================================================

GLOBALS_SECTION
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

DATA_SECTION
  // Read data file
  init_int nc
  init_matrix Cdata(1,2,1,nc)  // Year | C
  init_int ni
  init_int nIdx
  init_matrix Idata(1,3,1,ni)  // Year | I
  
  // Vectors
  ivector Cyear(1,nc)
  vector      C(1,nc)

  ivector Iyear(1,ni)
  ivector   Idx(1,ni)
  vector      I(1,ni)		
  ivector     X(1,ni)  // years with abundance index: 1995 | 1998 | ...
  vector   logI(1,ni)
  
  // Constants
  number halfnlog2pi
  

  number nReg

  // Switch to control file
  !! string run_name = string(adprogram_name);
  !! if(option_match(argc,argv,"-ind") > -1){
  !!   run_name = argv[option_match(argc,argv,"-ind") + 1];
  !!   run_name = run_name.substr(0, run_name.rfind("."));}
  
  // Read control file (phase, lower, upper, init)
  !! change_datafile_name((adstring)run_name.c_str() + ".ctl");

  init_vector _r_plui(1,4)
  init_vector _k_plui(1,4)
  init_vector _p_plui(1,4)
  init_vector _a_plui(1,4)

  init_ivector qPh(1,nIdx)
  init_vector  qLo(1,nIdx)
  init_vector  qHi(1,nIdx)
  init_vector  qPr(1,nIdx)
  
  init_ivector sPh(1,nIdx)
  init_vector  sLo(1,nIdx)
  init_vector  sHi(1,nIdx)
  init_vector  sPr(1,nIdx)
 
  // Switch to prior file
  !! change_datafile_name((adstring)run_name.c_str() + ".prr");

  // Read prior file (wt, mean, sd)
  init_vector r_prior(1,4)
  init_vector k_prior(1,4)
  init_vector p_prior(1,4)
  init_vector a_prior(1,4)
  init_vector q_prior(1,4)
  init_vector s_prior(1,4)
  
PARAMETER_SECTION
  // Estimated
  !! phz = (int) _r_plui[1];
  !! lb  =       _r_plui[2];
  !! ub  =       _r_plui[3];
  init_bounded_number _r(lb,ub,phz)
  !! phz = (int) _k_plui[1];
  !! lb  =       _k_plui[2];
  !! ub  =       _k_plui[3];
  init_bounded_number _k(lb,ub,phz)
  !! phz = (int) _a_plui[1];
  !! lb  =       _a_plui[2];
  !! ub  =       _a_plui[3];
  init_bounded_number _a(lb,ub,phz)
  !! phz = (int) _p_plui[1];
  !! lb  =       _p_plui[2];
  !! ub  =       _p_plui[3];
  init_bounded_number _p(lb,ub,phz)
  
  //change from log
  //init_bounded_number_vector logq(1,nIdx,qLo,qHi,qPh)  
  //init_bounded_number_vector logs(1,nIdx,sLo,sHi,sPh)  
  init_bounded_number_vector _q(1,nIdx,qLo,qHi,qPh)  
  init_bounded_number_vector _s(1,nIdx,sLo,sHi,sPh)  

  // Derived
  sdreport_number r
  sdreport_number k
  sdreport_number a
  sdreport_number p
  sdreport_vector q(1,nIdx)
  sdreport_vector s(1,nIdx)
  sdreport_number cnow
  sdreport_number bnow
  sdreport_number fnow
  sdreport_number msy
  sdreport_number bmsy
  sdreport_number fmsy
  sdreport_number cmsy
  sdreport_number bbmsy
  sdreport_number ffmsy
  sdreport_number bk
  sdreport_number fr
  sdreport_number bratio
  sdreport_number fratio
  sdreport_number slopeb;
  sdreport_number slopef;
  sdreport_vector ll(1,nIdx)
  
  number SS
  number s_
  vector nii(1,2)
  
  number _bratio
  number _fratio

  number  xy  
  number  x   
  number  y 
  number  xx

  // Updated
  vector B(1,nc+1)
  vector Bfit(1,ni)
  vector Ifit(1,ni)
  vector  RSS(1,nIdx)
    
  // Report
  matrix summary(1,nc,1,7)  // Year | B | C | I | Ifit | stockHat | stock.

  // likelihood profile
  likeprof_number lpr 

  // Objfun
  objective_function_value neglogL

PRELIMINARY_CALCS_SECTION
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
  lpr   =_r;
  
PROCEDURE_SECTION
  get_fit();
  get_neglogL();
  if(mceval_phase())
    write_mcmc();

  // likelihood profile 
  lpr=_r;
  
REPORT_SECTION
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
 
FUNCTION get_fit
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
   

FUNCTION get_neglogL
 
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
  if (r_prior[1]>1) neglogL -= dnorm(r, r_prior[2], r_prior[3])/dnorm(r_prior[2], r_prior[2], r_prior[3]);
  if (k_prior[1]>1) neglogL -= dnorm(k, k_prior[2], k_prior[3])/dnorm(k_prior[2], k_prior[2], k_prior[3]);
  if (p_prior[1]>1) neglogL -= dnorm(p, p_prior[2], p_prior[3])/dnorm(p_prior[2], p_prior[2], p_prior[3]);
  if (a_prior[1]>1) neglogL -= dnorm(a, a_prior[2], a_prior[3])/dnorm(a_prior[2], a_prior[2], a_prior[3]);
  //for (i=1; i<=nIdx; i++){
  //  neglogL += q_prior[i]*dnorm(q[i],    q_prior[i,2],     q_prior[i,3]);
  //  neglogL += s_prior[i]*dnorm(s[i],    s_prior[i,2],     s_prior[i,3]);}
   

FUNCTION get_summary
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

FUNCTION write_mcmc
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

FUNCTION write_priors
  priors<<r_prior<<endl;
	
FUNCTION write_bounds
  bounds<<q	
        <<s<<endl;	

FUNCTION write_ll
   for (int j=1; j<=nIdx; j++) 
       ll[j]=0.0; 

    for (int j=1; j<=ni; j++){
       s_    = _s[Idx[j]];
       ll[Idx[j]] += log(s_) +pow(log(I[j])-log(Ifit[j]),2.0)/(2*s_*s_);
       }

  lls<<"# ll"<<endl<<ll<<endl;

TOP_OF_MAIN_SECTION
  arrmblsize = 40000000L;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(3000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(200000);
  gradient_structure::set_MAX_NVAR_OFFSET(10000);

