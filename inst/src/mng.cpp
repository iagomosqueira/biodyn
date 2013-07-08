SEXP mng(SEXP stock, SEXP harvest, SEXP refpt, SEXP rtn, int nStep=3, int nReg=3){
  
  FLQuant stk(stock);
  FLQuant hvt(harvest);
  FLQuant ref(refpt);
  FLQuant rtn(rtn);
 
  int i, j, k, l, m;
  for(i=1; i<=stk.niters();     l++)
    for (j=1; j<=stk.nareas();    j++)
      for (k=1; k<=syk.nseasons();  k++) 
        for (l=1; l<=stk.nunits();    l++){
		  
		  // alsolute values
		  // cnow 
		  rtn[1,1,l,k,j,i]=C[nc];
		  // fnow 
		  rtn[1,2,l,k,j,i]=C[nc]/B[nc];
		  // bnow 
		  rtn[1,3,l,k,j,i]=B[nc];
			    
		  // relative values
		  //cmsy 
		  rtn[1,1,l,k,j,i]=C[nc]/msy;
		  //bbmsy
		  rtn[1,1,l,k,j,i]=bnow/bmsy;
		  //ffmsy
		  rtn[1,1,l,k,j,i]=fnow/fmsy;
		  //bk   
		  rtn[1,1,l,k,j,i]=bnow/k;
		  //fr   
		  rtn[1,1,l,k,j,i]=fnow/r;
			  
		  // "last time"
		  //bratio =0.0;
		  //fratio =0.0;
		  //_bratio=0.0;
		  //_fratio=0.0;
			  
		  for (int i=nc; i>nc-3; i--){
			//bratio+=B[i];
			//fratio+=C[i]/B[i];
			
			//_bratio+=B[i-3];
			//_fratio+=C[i-3]/B[i-3];
			}
		  //bratio=bratio/_bratio;
		  //fratio=fratio/_fratio;

		  // slopes
		  xy=x=y=xx=0.0; 
		  for (int i=nc; i>nc-nReg; i--){
			x +=i;
			xx+=i*i;  
			y +=B[i];  
			xy+=i*B[i];  
			}
		  slopeb = (nReg*xy - x*y)/(nReg*xx - x*2.0);
			   
		  xy=0.0; 
		  x =0.0;
		  y =0.0; 
		  xx=0.0; 
		  for (int i=nc; i>nc-nReg; i--){
			x +=i;
			xx+=i*i;  
			y +=C[i]/B[i];  
			xy+=i*C[i]/B[i];  
			}
		  slopef = (nReg*xy - x*y)/(nReg*xx - x*2.0);
		  }
			   
  return rtn;
  }
