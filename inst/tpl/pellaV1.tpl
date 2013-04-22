//=========================================================================================================================
// File:        pella.tpl
// Model:       Pella-Tomlinson model, with Binit=k*a
// Parameters:  r, k, a, p, q, sigma
// Fitted data: Abundance index
// Likelihood:  Log-transformed normal
// References:  Polacheck et al. (1993)
// Notes:       q and sigma are free parameters, to allow full uncertainty in MCMC analysis
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

DATA_SECTION
  // Read data file
  init_int nc
  init_matrix Cdata(1,2,1,nc)  // Year | C
  init_int ni
  init_matrix Idata(1,2,1,ni)  // Year | I
  
  // Vectors
  ivector Cyear(1,nc)
  ivector Iyear(1,ni)
  vector  C(1,nc)
  vector  I(1,ni)		
  ivector X(1,ni)  // years with abundance index: 1995 | 1998 | ...
  vector logI(1,ni)
  
  // Constants
  number halfnlog2pi
  
  //checks
  number q1
  number q2
  number sigma2

  // Switch to control file
  !! string run_name = string(adprogram_name);
  !! if(option_match(argc,argv,"-ind") > -1){
  !!   run_name = argv[option_match(argc,argv,"-ind") + 1];
  !!   run_name = run_name.substr(0, run_name.rfind("."));}
  
  // Read control file (phase, lower, upper, init)
  !! change_datafile_name((adstring)run_name.c_str() + ".ctl");

  init_vector logr_plui(1,4)
  init_vector logk_plui(1,4)
  init_vector logp_plui(1,4)
  init_vector loga_plui(1,4)
  init_vector logq_plui(1,4)
  init_vector logsigma_plui(1,4)

  // Switch to prior file
  !! change_datafile_name((adstring)run_name.c_str() + ".prr");

  // Read prior file (wt, mean, sd)
  init_vector r_prior(1,4)
  init_vector k_prior(1,4)
  init_vector p_prior(1,4)
  init_vector a_prior(1,4)
  init_vector q_prior(1,4)
  init_vector sigma_prior(1,4)
  
PARAMETER_SECTION
  // Estimated
  !! phz = (int) logr_plui[1];
  !! lb  =       logr_plui[2];
  !! ub  =       logr_plui[3];
  init_bounded_number logr(lb,ub,phz)
  !! phz = (int) logk_plui[1];
  !! lb  =       logk_plui[2];
  !! ub  =       logk_plui[3];
  init_bounded_number logk(lb,ub,phz)
  !! phz = (int) loga_plui[1];
  !! lb  =       loga_plui[2];
  !! ub  =       loga_plui[3];
  init_bounded_number loga(lb,ub,phz)
  !! phz = (int) logp_plui[1];
  !! lb  =       logp_plui[2];
  !! ub  =       logp_plui[3];
  init_bounded_number logp(lb,ub,phz)
  !! phz = (int) logq_plui[1];
  !! lb  =       logq_plui[2];
  !! ub  =       logq_plui[3];
  init_bounded_number logq(lb,ub,phz)
  !! phz = (int) logsigma_plui[1];
  !! lb  =       logsigma_plui[2];
  !! ub  =       logsigma_plui[3];
  init_bounded_number logsigma(lb,ub,phz)
  
  // Derived
  sdreport_number r
  sdreport_number k
  sdreport_number a
  sdreport_number p
  sdreport_number q
  sdreport_number B2
  sdreport_number BBMSY
  sdreport_number HHMSY
  number sigma
  number SS
  
  // Updated
  vector B(1,nc)
  vector Bfit(1,ni)
  vector Ifit(1,ni)
  number RSS
  
  // Report
  matrix summary(1,nc,1,5)  // Year | B | C | I | Ifit
  
  // Objfun
  objective_function_value neglogL

PRELIMINARY_CALCS_SECTION
  // Data
  Cyear = (ivector) row(Cdata,1);
  C = row(Cdata,2);
  Iyear = (ivector) row(Idata,1);
  I = row(Idata,2);
  X = Iyear - Cyear[1] + 1;
  halfnlog2pi = 0.5*ni*log(2*pi);
  logI = log(I);

  // Parameters
  logr = logr_plui[4];
  logk = logk_plui[4];
  loga = loga_plui[4];
  logp = logp_plui[4];
  logq = logq_plui[4];
  logsigma = logsigma_plui[4];

PROCEDURE_SECTION
  get_fit();
  get_neglogL();
  if(mceval_phase())
    write_mcmc();
  
  //write_priors();

REPORT_SECTION
  summary.initialize();
  get_summary();
 
  report<<setprecision(12)
        <<"# r"      <<endl<<r      <<endl
        <<"# k"      <<endl<<k      <<endl
        <<"# b0"     <<endl<<a      <<endl
        <<"# p"      <<endl<<p      <<endl
        <<"# q"      <<endl<<q      <<endl
        <<"# sigma"  <<endl<<sigma  <<endl
        <<"# RSS"    <<endl<<RSS    <<endl
        <<"# neglogL"<<endl<<neglogL<<endl<<endl;
  report<<setprecision(12)
        <<"# Model summary"<<endl
        <<" Year Biomass Catch Index IndexFit"<<endl
        <<summary<<endl;

FUNCTION get_fit
  r = mfexp(logr);
  k = mfexp(logk);
  a = mfexp(loga);
  p = mfexp(logp);

  q = mfexp(logq);
  
  //Check
  //q1=q2=0.0;
  //for(int i=1; i<=ni; i++) {
  //  int j  =X[i];
  //  q1 +=B(j)*I(i);
  //  q2 +=B[j]*B[j];}
  //q2=q1/q2;

  B[1] = k * a;
  for(int t=1; t<=nc-1; t++)
    B[t+1] = sfabs(B[t] + r/p*B[t]*(1-pow(B[t]/k,p)) - C[t]);
  Ifit = q*B(X);

  sigma = mfexp(logsigma);
  //SS    = 0.0;
  //for(int i=1; i<=ni; i++) 
  //   SS   +=pow(log(I(i))-log(Ifit(i)),2);
  //sigma2=log(pow((SS/ni),.5));

  B2   =B[nc-1];
  BBMSY=B2/(k*pow((1/(1+p)),(1/p)));
  HHMSY=C[nc-1]/B[nc-1]/(r*k*pow((1/(1+p)),(1/p+1))/ k*pow((1/(1+p)),(1/p)));

FUNCTION get_neglogL
  RSS     = norm2(logI-log(Ifit));
  neglogL = halfnlog2pi + ni*log(sigma) + RSS/(2*sigma*sigma);

  // weighted likelihood priors
  if (r_prior[1]>0) 
      neglogL += log(r_prior[1]*dnorm(r,    r_prior[2],     r_prior[3]));
  if (k_prior[1]>0)
      neglogL += log(k_prior[1]*dnorm(k,    k_prior[2],     k_prior[3]));
  if (p_prior[1]>0)
      neglogL += log(p_prior[1]*dnorm(p,    p_prior[2],     p_prior[3]));
  if (a_prior[1]>0) 
      neglogL += log(a_prior[1]*dnorm(a,    a_prior[2],     a_prior[3]));
  if (q_prior[1]>0)
      neglogL += log(q_prior[1]*dnorm(q,    q_prior[2],     q_prior[3]));
  if (sigma_prior[1]>0)
      neglogL += log(sigma_prior[1]*dnorm(sigma,sigma_prior[2], sigma_prior[3]));
   


FUNCTION get_summary
  summary.colfill(1,(dvector)Cyear);
  summary.colfill(2,B);
  summary.colfill(3,C);
  for(int i=1; i<=ni; i++)  // allow missing years in abundance index
  {
    summary(X[i],4) = I[i];
    summary(X[i],5) = Ifit[i];
  }

FUNCTION write_mcmc
  mcmc_iteration++;
  // Parameters
  if(mcmc_iteration == 1){
    mcmc_par<<"neglogL,r,k,b0,p,q,sigma"<<endl;
    mcmc_par<<neglogL<<","
          <<r      <<","
          <<k      <<","
          <<a      <<","
          <<p      <<","
          <<q      <<","
          <<sigma  <<endl;}
  // Biomass
  if(mcmc_iteration == 1)
  {
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
  priors<<p_prior<<endl;	
