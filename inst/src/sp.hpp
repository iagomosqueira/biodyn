#ifndef _INC_SP
#define _INC_SP

#include <FLCoreClasses.hpp>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <FLCoreClasses.hpp>

#define const_nan 0.0

typedef enum tagFLRConstSP 
   {
   FLRConst_fletcher =1,
   FLRConst_fox      =2,
   FLRConst_genfit   =3,
   FLRConst_gulland  =4,
   FLRConst_logistic =5,
   FLRConst_pellat   =6,
   FLRConst_schaefer =7,
   FLRConst_shepherd =8  
   } FLRConstSP;

FLRConstSP getSP(int        i);
int        getSP(FLRConstSP i);

SEXP msy( int model, SEXP params);
SEXP fmsy(int model, SEXP params);
SEXP bmsy(int model, SEXP params);
SEXP sp(  int model, SEXP params);

class SP
{
public:        
   SP(int);
   SP(int,FLPar);
   
   FLQuant  msy(int, FLPar);
   FLQuant fmsy(int, FLPar);
   FLQuant bmsy(int, FLPar);
   FLQuant sp(  int, FLPar);

   SEXP  msy(void);
   SEXP fmsy(void);
   SEXP bmsy(void);
   SEXP func(SEXP);


  ~SP(void);      
 
protected: 
   FLRConstSP  model;
   FLPar       par;
  
   void unalloc(void);
   };                  

#endif /* _INC_SP */




