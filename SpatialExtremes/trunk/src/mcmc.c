#include "header.h"

SEXP gibbs(SEXP n, SEXP np, SEXP thin, SEXP init,
	   SEXP propsd, SEXP f, SEXP rho){

  int i,j,k,nr;
  int nn = INTEGER(n)[0], nnp = INTEGER(np)[0], thinn = INTEGER(thin)[0];
  double prop, prop_ratio, acc_prob, post_ratio;
  double *crow, *prow;
  SEXP ans, nacc, nex, mc, current, dpst_lower, dpst_upper;

  nr = 1 + ftrunc(nn/thinn);
  crow = (double *)R_alloc(nnp, sizeof(double));
  prow = (double *)R_alloc(nnp, sizeof(double));
  PROTECT(current = allocVector(REALSXP, nnp));
  PROTECT(nacc = allocVector(INTSXP, nnp));
  PROTECT(nex = allocVector(INTSXP, nnp));
  PROTECT(mc = allocVector(REALSXP, nr * nnp));
  PROTECT(ans = allocVector(VECSXP, 3));
  PROTECT(dpst_lower = allocVector(REALSXP, 1));
  PROTECT(dpst_upper = allocVector(REALSXP, 1));

  for(i=0;i<nnp;i++) {
    prow[i] = REAL(init)[i];
    REAL(mc)[i] = REAL(init)[i];
    INTEGER(nex)[i] = INTEGER(nacc)[i] = 0;
  }

  RANDIN;
  for(i=0;i<nn;i++) {
    for(j=0;j<nnp;j++) {     

      if (j < 2){
	prop = rlnorm(log(prow[j]), REAL(propsd)[j]);
	prop_ratio = prop / prow[j];
      }

      else{
	prop = rnorm(prow[j], REAL(propsd)[j]);
	prop_ratio = 1;
      }
      
      for(k=0;k<nnp;k++) {
     	if (k < j) 
	  REAL(current)[k] = crow[k];

        else
	  REAL(current)[k] = prow[k];
      }

      defineVar(install("x"), current, rho);
      dpst_lower = eval(f, rho);
      REAL(current)[j] = prop;

      defineVar(install("x"), current, rho);
      dpst_upper = eval(f, rho);

      post_ratio = exp(REAL(dpst_upper)[0] - REAL(dpst_lower)[0]);
      
      if(!R_FINITE(REAL(dpst_upper)[0]))
        INTEGER(nex)[j] = INTEGER(nex)[j] + 1;

      acc_prob = fmin2(1.0, prop_ratio * post_ratio);

      if(!R_FINITE(acc_prob)) {
        acc_prob = 0;
        warning("NaN returned for posterior density");
      }

      if (runif(0, 1) < acc_prob) {
        crow[j] = prop;
	INTEGER(nacc)[j] = INTEGER(nacc)[j] + 1;
      }

      else crow[j] = prow[j];
    }

    if( ((i+1) % thinn) == 0)
      for(j=0;j<nnp;j++) {
	REAL(mc)[(i+1)/thinn * nnp + j] = crow[j];
      }
    
    for (j=0;j<nnp;j++)
      prow[j] = crow[j];
  }

  RANDOUT;
  SET_VECTOR_ELT(ans, 0, mc);
  SET_VECTOR_ELT(ans, 1, nacc);
  SET_VECTOR_ELT(ans, 2, nex);
  UNPROTECT(7);
  return(ans);
}
