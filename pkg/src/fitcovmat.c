#include "header.h"

void fitcovmat2d(double *cov11, double *cov12, double *cov22,
		 int *nPairs, double *distVec, double *extcoeff,
		 double *weights, double *ans){
  
  int i;
  double *mahalDist, nonfeas = 0.0;

  mahalDist = (double *)R_alloc(*nPairs, sizeof(double));

  nonfeas = mahalDistFct(distVec, *nPairs, cov11, cov12, cov22,
			 mahalDist);

  for (i=0;i<*nPairs;i++)
    *ans += R_pow_di((2 * pnorm(mahalDist[i] / 2, 0.0, 1.0, 1, 0)
		      - extcoeff[i]) / weights[i], 2);

  if (nonfeas != 0.0)
    *ans = nonfeas * *ans;
  return;
}


void fitcovmat3d(double *cov11, double *cov12, double *cov13,
		 double *cov22, double *cov23, double *cov33,
		 int *nPairs, double *distVec, double *extcoeff,
		 double *weights, double *ans){
  
  int i;
  double *mahalDist, nonfeas = 0.0;

  mahalDist = (double *)R_alloc(*nPairs, sizeof(double));

  nonfeas = mahalDistFct3d(distVec, *nPairs, cov11, cov12, cov13, 
			   cov22, cov23, cov33, mahalDist);

  for (i=0;i<*nPairs;i++)
    *ans += R_pow_di((2 * pnorm(mahalDist[i] / 2, 0.0, 1.0, 1, 0) -
		      extcoeff[i]) / weights[i], 2);

  if (nonfeas != 0.0)
    *ans = nonfeas * *ans;
  
  return;
}

void fitcovariance(int *covmod, double *sill, double *range, double *smooth,
		   int *nPairs, double *dist, double *extcoeff,
		   double *weights, double *ans){
  
  int i;
  double *rho, nonfeas = 0.0;

  rho = (double *)R_alloc(*nPairs, sizeof(double));

  switch (*covmod){
  case 1:
    nonfeas = whittleMatern(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  case 2:
    nonfeas = cauchy(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  case 3:
    nonfeas = powerExp(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  }

  for (i=0;i<*nPairs;i++)
    *ans += R_pow_di((1 + sqrt(1 - 0.5 * (rho[i] + 1)) -
		      extcoeff[i]) / weights[i], 2);

  if (nonfeas != 0.0)
    *ans = nonfeas * *ans;

  return;
}
