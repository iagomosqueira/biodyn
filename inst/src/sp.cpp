#include <sp.hpp>
#include <Rcpp.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <FLCoreClasses.hpp>

FLRConstSP getSP(int i)
   {
   switch(i) {
   	case 1:  return FLRConst_fletcher;
   	case 2:  return FLRConst_fox;
   	case 3:  return FLRConst_genfit;
   	case 4:  return FLRConst_gulland;
   	case 5:  return FLRConst_logistic;
   	case 6:  return FLRConst_pellat;
   	case 7:  return FLRConst_schaefer;
   	case 8:  return FLRConst_shepherd;  
      }

   return FLRConst_pellat;}

int getSP(FLRConstSP i)
   {
   switch(i) {
	case FLRConst_fletcher: return 1;
   	case FLRConst_fox:    	return 2;
   	case FLRConst_genfit:   return 3;
   	case FLRConst_gulland:  return 4;
   	case FLRConst_logistic: return 5;
   	case FLRConst_pellat:   return 6;
   	case FLRConst_schaefer: return 7;
   	case FLRConst_shepherd: return 8;
        }

   return 6;}


 // MSY
double msyFox(double r, double k) { return r*(k*exp(-1))*(1-(log(k)-1)/log(k));}
double msySchaefer(double r,double k) {return r*k/4;}
double msyLogistic(double msy,double k){return msy;}
double msyPellaT(double r,double k, double p){return r*k*(1/(1+p))^(1/p+1);}
double msyGenfit(double r, double k, double p){ return r*k*(1/(1+p))^(1/p+1);}
double msyGulland (double r, double k){return (r*k^2)/4;}
double msyFletcher(double msy, double k, double p){return msy;}
double msyShepherd(double r, double k, double m) {
  double aPrime=r/m - 1.0;
  double Bmax  =k*aPrime;
  double _bmsy = 0; //bmsy("shepherd" double p,aram);
  
  return aPrime*m*_bmsy*(1-_bmsy/Bmax)/(1+aPrime)^.5;}

// BMSY
double bmsyFox (double r, double k){return k*exp(-1);}
double bmsySchaefer(double r, double k){return k/2;}
double bmsyLogistic(double r, double k){return k/2;}
double bmsyPellaT(double r, double k, double p){return k*(1/(1+p))^(1/p);}
double bmsyGenfit(double r, double k, double p){return k*(1/(1+p))^(1/p);}
double bmsyGulland(double r, double k){return k/2;}
double bmsyFletcher(double msy, double k, double p) {return k*(1/(p+1)^(1/(p)));}
double bmsyShepherd(double r, double k, double m) {
  double aPrime = r/m - 1;
  double Bmax  = k*aPrime;
  
  return Bmax*((1+aPrime)^.5-1)/aPrime;}

// FMSY
double fmsyPellaT(double r, double k, double p)   {return r*(1/(1+p));}
double fmsyGenfit(double r, double k, double p)   {return r*(1/(1+p));}
double fmsyFox(double r, double k)       {return r*(1-(log(k)-1)/log(k));}
double fmsySchaefer(double r, double k) {return r/2;}
double fmsyShepherd(double r, double k, double m) {return msyShepherd(r,k,m)/bmsyShepherd(r,k,m);}
double fmsyGulland(double r, double k)  {return r*k/2;}
double fmsyFletcher(double msy, double k, double p) {return msyFletcher(msy,k,p)/bmsyFletcher(msy,k,p);}
double fmsyLogistic(double msy, double k) {
	double r=4*msy/k;
	return r/2;}

//Production Functions
double spFox(double biomass, double r, double k){
	return r*biomass*(1-log(biomass)/log(k));}
double spSchaefer(double biomass, double r, double k){
     return r*biomass*(1-biomass/k);}
double spPellaT(double biomass, double r, double k, double p){
    double a=biomass*r/p;
    double b=biomass/k;
    double c=b**p;
    return  a*(1-c);}
double spShepherd(double biomass, double r, double k, double m){
    return r*biomass/(1+biomass/k)-m*biomass;}
double spGulland(double biomass, double r, double k){
    return r*biomass*(k-biomass);}
double spFletcher(double biomass, double msy, double k, double p) {
    p      = p+1;
    double lambda = (p**(p/(p-1)))/(p-1);
    return lambda*msy*(biomass/k)-lambda*msy*(biomass/k)**p;}
double spLogistic(double biomass, double msy, double k){
    double r=4*msy/k;
    return r*biomass*(1-biomass/k);}
double spGenfit(double biomass, double r, double k, double p){
    return r*biomass*(1-biomass/k);}

  
FLQuant msy(int model, SEXP params){
   FLQuant par(params);
   FLQuant rtn(1,1,par.minyr(),par.maxyr(),par.nareas(),par.nseasons(),par.nunits(),par.niters());

   int i, j, k, l, m; 
   for(m=1; m <= par.niters(); m++){
     for (l = 1; l <= par.nareas(); l++){
       for (k = 1; k <= par.nseasons(); k++) {
	     for (j = 1; j <= par.nunits(); j++){
	       for (i = par.minyr(); i <= par.maxyr(); i++){

          switch(getSP(model)) {
				case FLRConst_fletcher: msyFletcher(par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_fox:      msyFox(     par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_genfit:   msyGenfit(  par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_gulland:  msyGulland( par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_logistic: msyLogistic(par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_pellat:   msyPellaT(  par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_schaefer: msySchaefer(par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_shepherd: msyShepherd(par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			    }}}}}}

   return rtn;}

FLQuant bmsy(int model, SEXP params){
	   FLQuant par(params);
	   FLQuant rtn(1,1,par.minyr(),par.maxyr(),par.nareas(),par.nseasons(),par.nunits(),par.niters());

	   int i, j, k, l, m;
	   for(m=1; m <= par.niters(); m++){
	     for (l = 1; l <= par.nareas(); l++){
	       for (k = 1; k <= par.nseasons(); k++) {
		     for (j = 1; j <= par.nunits(); j++){
		       for (i = par.minyr(); i <= par.maxyr(); i++){

   			  switch(getSP(model)) {
				case FLRConst_fletcher: rtn(1,i,j,k,l,m) = bmsyFletcher(par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_fox:      rtn(1,i,j,k,l,m) = bmsyFox(     par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_genfit:   rtn(1,i,j,k,l,m) = bmsyGenfit(  par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_gulland:  rtn(1,i,j,k,l,m) = bmsyGulland( par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_logistic: rtn(1,i,j,k,l,m) = bmsyLogistic(par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_pellat:   rtn(1,i,j,k,l,m) = bmsyPellaT(  par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_schaefer: rtn(1,i,j,k,l,m) = bmsySchaefer(par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_shepherd: rtn(1,i,j,k,l,m) = bmsyShepherd(par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
                }}}}}}

   return rtn;}

FLQuant fmsy(int model, SEXP params){
	   FLQuant par(params);
	   FLQuant rtn(1,1,par.minyr(),par.maxyr(),par.nareas(),par.nseasons(),par.nunits(),par.niters());

	   int i, j, k, l, m;
	   for(m=1; m <= par.niters(); m++){
	     for (l = 1; l <= par.nareas(); l++){
	       for (k = 1; k <= par.nseasons(); k++) {
		     for (j = 1; j <= par.nunits(); j++){
		       for (i = par.minyr(); i <= par.maxyr(); i++){
			       			   
   			   switch(getSP(model)) {
				case FLRConst_fletcher: rtn(1,i,j,k,l,m) = fmsyFletcher(par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_fox:      rtn(1,i,j,k,l,m) = fmsyFox(     par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_genfit:   rtn(1,i,j,k,l,m) = fmsyGenfit(  par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_gulland:  rtn(1,i,j,k,l,m) = fmsyGulland( par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_logistic: rtn(1,i,j,k,l,m) = fmsyLogistic(par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_pellat:   rtn(1,i,j,k,l,m) = fmsyPellaT(  par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_schaefer: rtn(1,i,j,k,l,m) = fmsySchaefer(par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_shepherd: rtn(1,i,j,k,l,m) = fmsyShepherd(par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
                }}}}}}

   return rtn;}

FLQuant sp(int model, SEXP params, SEXP biomass)
   {
	FLQuant bio(biomass);
	FLQuant par(params);
	FLQuant rtn(1,1,par.minyr(),par.maxyr(),par.nareas(),par.nseasons(),par.nunits(),par.niters());

   int i, j, k, l, m; 
   int i, j, k, l, m;
   for(m=1; m <= bio.niters(); m++){
     for (l = 1; l <= bio.nareas(); l++){
       for (k = 1; k <= bio.nseasons(); k++) {
	     for (j = 1; j <= bio.nunits(); j++){
	       for (i = bio.minyr(); i <= bio.maxyr(); i++){
			       			   
	   		   switch(getSP(model)) {
	   			case FLRConst_fletcher: rtn(1,i,j,k,l,m) = spFletcher(bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_fox:      rtn(1,i,j,k,l,m) = spFox(     bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_genfit:   rtn(1,i,j,k,l,m) = spGenfit(  bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_gulland:  rtn(1,i,j,k,l,m) = spGulland( bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_logistic: rtn(1,i,j,k,l,m) = spLogistic(bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_pellat:   rtn(1,i,j,k,l,m) = spPellaT(  bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_schaefer: rtn(1,i,j,k,l,m) = spSchaefer(bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_shepherd: rtn(1,i,j,k,l,m) = spShepherd(bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
                }}}}}}

   return rtn;}


SEXP SP::msy(void){
   FLQuant rtn(1,1,par.minyr(),par.maxyr(),par.nareas(),par.nseasons(),par.nunits(),par.niters());

   int i, j, k, l, m;
   for(m=1; m <= par.niters(); m++){
     for (l = 1; l <= par.nareas(); l++){
       for (k = 1; k <= par.nseasons(); k++) {
	     for (j = 1; j <= par.nunits(); j++){
	       for (i = par.minyr(); i <= par.maxyr(); i++){

          switch(getSP(model)) {
				case FLRConst_fletcher: msyFletcher(par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_fox:      msyFox(     par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_genfit:   msyGenfit(  par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_gulland:  msyGulland( par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_logistic: msyLogistic(par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_pellat:   msyPellaT(  par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_schaefer: msySchaefer(par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_shepherd: msyShepherd(par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			    }}}}}}

   return rtn;}

SEXP SP::bmsy(void){
	   FLQuant rtn(1,1,par.minyr(),par.maxyr(),par.nareas(),par.nseasons(),par.nunits(),par.niters());

	   int i, j, k, l, m;
	   for(m=1; m <= par.niters(); m++){
	     for (l = 1; l <= par.nareas(); l++){
	       for (k = 1; k <= par.nseasons(); k++) {
		     for (j = 1; j <= par.nunits(); j++){
		       for (i = par.minyr(); i <= par.maxyr(); i++){

   			  switch(getSP(model)) {
				case FLRConst_fletcher: rtn(1,i,j,k,l,m) = bmsyFletcher(par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_fox:      rtn(1,i,j,k,l,m) = bmsyFox(     par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_genfit:   rtn(1,i,j,k,l,m) = bmsyGenfit(  par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_gulland:  rtn(1,i,j,k,l,m) = bmsyGulland( par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_logistic: rtn(1,i,j,k,l,m) = bmsyLogistic(par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_pellat:   rtn(1,i,j,k,l,m) = bmsyPellaT(  par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_schaefer: rtn(1,i,j,k,l,m) = bmsySchaefer(par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_shepherd: rtn(1,i,j,k,l,m) = bmsyShepherd(par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
                }}}}}}

   return rtn;}

SEXP SP::fmsy(void){
	   FLQuant rtn(1,1,par.minyr(),par.maxyr(),par.nareas(),par.nseasons(),par.nunits(),par.niters());

	   int i, j, k, l, m;
	   for(m=1; m <= par.niters(); m++){
	     for (l = 1; l <= par.nareas(); l++){
	       for (k = 1; k <= par.nseasons(); k++) {
		     for (j = 1; j <= par.nunits(); j++){
		       for (i = par.minyr(); i <= par.maxyr(); i++){

   			   switch(getSP(model)) {
				case FLRConst_fletcher: rtn(1,i,j,k,l,m) = fmsyFletcher(par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_fox:      rtn(1,i,j,k,l,m) = fmsyFox(     par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_genfit:   rtn(1,i,j,k,l,m) = fmsyGenfit(  par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_gulland:  rtn(1,i,j,k,l,m) = fmsyGulland( par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_logistic: rtn(1,i,j,k,l,m) = fmsyLogistic(par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_pellat:   rtn(1,i,j,k,l,m) = fmsyPellaT(  par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_schaefer: rtn(1,i,j,k,l,m) = fmsySchaefer(par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_shepherd: rtn(1,i,j,k,l,m) = fmsyShepherd(par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
                }}}}}}

   return rtn;}

SEXP SP::func(SEXP biomass)
   {
	FLQuant bio(biomass);
	FLQuant rtn(1,1,par.minyr(),par.maxyr(),par.nareas(),par.nseasons(),par.nunits(),par.niters());

   int i, j, k, l, m;
   int i, j, k, l, m;
   for(m=1; m <= bio.niters(); m++){
     for (l = 1; l <= bio.nareas(); l++){
       for (k = 1; k <= bio.nseasons(); k++) {
	     for (j = 1; j <= bio.nunits(); j++){
	       for (i = bio.minyr(); i <= bio.maxyr(); i++){

	   		   switch(getSP(model)) {
	   			case FLRConst_fletcher: rtn(1,i,j,k,l,m) = spFletcher(bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_fox:      rtn(1,i,j,k,l,m) = spFox(     bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_genfit:   rtn(1,i,j,k,l,m) = spGenfit(  bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_gulland:  rtn(1,i,j,k,l,m) = spGulland( bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_logistic: rtn(1,i,j,k,l,m) = spLogistic(bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_pellat:   rtn(1,i,j,k,l,m) = spPellaT(  bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
			   	case FLRConst_schaefer: rtn(1,i,j,k,l,m) = spSchaefer(bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m));                  break;
			   	case FLRConst_shepherd: rtn(1,i,j,k,l,m) = spShepherd(bio(1,i,j,k,l,m),par(1,i,j,k,l,m),par(2,i,j,k,l,m),par(3,i,j,k,l,m)); break;
                }}}}}}

   return rtn;}

