#include "header.h"

/* These functions estimate the spatial dependence parameters by using
   a least square optimization. More precisely, it compares the
   expected extremal coefficients to the "observed" ones.
*/

void fitcovmat2d(double *cov11, double *cov12, double *cov22,
		 int *nPairs, double *distVec, double *extcoeff,
		 double *weights, double *ans){
  
  int i;
  double *mahalDist, res;

  mahalDist = (double *)R_alloc(*nPairs, sizeof(double));

  *ans = - mahalDistFct(distVec, *nPairs, cov11, cov12, cov22,
			mahalDist);

  if (*ans != 0.0){
    *ans = 1e50;
    return;
  }

  for (i=0;i<*nPairs;i++){
    res = 2 * pnorm(0.5 * mahalDist[i], 0, 1, 1, 0) - extcoeff[i];
    *ans +=  res * res / (weights[i] * weights[i]);
  }

  return;
}


void fitcovmat3d(double *cov11, double *cov12, double *cov13,
		 double *cov22, double *cov23, double *cov33,
		 int *nPairs, double *distVec, double *extcoeff,
		 double *weights, double *ans){
  
  int i;
  double *mahalDist, res;

  mahalDist = (double *)R_alloc(*nPairs, sizeof(double));

  *ans = - mahalDistFct3d(distVec, *nPairs, cov11, cov12, cov13, 
			  cov22, cov23, cov33, mahalDist);

  if (*ans != 0.0)
    return;

  for (i=0;i<*nPairs;i++){
    res = 2 * pnorm(0.5 * mahalDist[i], 0, 1, 1, 0) - extcoeff[i];
    *ans += res * res / (weights[i] * weights[i]);
  }

  return;
}

void fitcovariance(int *covmod, double *sill, double *range, double *smooth,
		   double *smooth2, int *nPairs, int *dim, double *dist,
		   double *extcoeff, double *weights, double *ans){
  
  int i;
  double *rho, res;

  rho = (double *)R_alloc(*nPairs, sizeof(double));

  switch (*covmod){
  case 1:
    *ans = -whittleMatern(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  case 2:
    *ans = -cauchy(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  case 3:
    *ans = -powerExp(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  case 4:
    *ans = -bessel(dist, *nPairs, *dim, *sill, *range, *smooth, rho);
    break;
  case 5:
    *ans = -caugen(dist, *nPairs, *sill, *range, *smooth, *smooth2, rho);
    break;
  }
  
  if (*ans != 0.0)
    return;

  for (i=0;i<*nPairs;i++){
    res = 1 + sqrt(0.5 - 0.5 * rho[i]) - extcoeff[i];
    *ans += res * res / (weights[i] * weights[i]);
  }

  return;
}

void fiticovariance(int *covmod, double *alpha, double *sill, double *range,
		    double *smooth, double *smooth2, int *nPairs, int *dim,
		    double *dist, double *extcoeff, double *weights, double *ans){
  /* This computes the least squares for the independent Schlather model */
  
  int i;
  double *rho, res;

  rho = (double *)R_alloc(*nPairs, sizeof(double));

  if (*alpha > 1){
    *ans = -*alpha * *alpha * MINF;
    return;
  }

  if (*alpha < 0){
    *ans = - (1 - *alpha) * (1 - *alpha) * MINF;
    return;
  }

  switch (*covmod){
  case 1:
    *ans = -whittleMatern(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  case 2:
    *ans = -cauchy(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  case 3:
    *ans = -powerExp(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  case 4:
    *ans = -bessel(dist, *nPairs, *dim, *sill, *range, *smooth, rho);
    break;
  case 5:
    *ans = -caugen(dist, *nPairs, *sill, *range, *smooth, *smooth2, rho);
    break;
  }
  
  if (*ans != 0.0)
    return;

  for (i=0;i<*nPairs;i++){
    res = 2 * *alpha + (1 - *alpha) * (1 + sqrt(0.5 - 0.5 * rho[i])) - extcoeff[i];
    *ans +=  res * res / (weights[i] * weights[i]);
  }

  return;
}

void fitgcovariance(int *covmod, double *sigma2, double *sigma2Bound, double *sill,
		    double *range, double *smooth, double *smooth2, int *nPairs,
		    int *dim, double *dist, double *extcoeff, double *weights,
		    double *ans){
  /* This computes the least squares for the geometric Gaussian model */
  
  int i;
  double *rho, res;

  rho = (double *)R_alloc(*nPairs, sizeof(double));

  *ans = -geomCovariance(dist, *nPairs, *dim, *covmod, *sigma2, *sigma2Bound,
			 *sill, *range, *smooth, *smooth2, rho);
    
  if (*ans != 0.0)
    return;

  for (i=0;i<*nPairs;i++){
    res = 2 * pnorm(0.5 * rho[i], 0.0, 1.0, 1, 0) - extcoeff[i];
    *ans += res * res / (weights[i] * weights[i]);
  }

  return;
}

void fitbrcovariance(double *range, double *smooth, int *nPairs,
		     double *dist, double *extcoeff, double *weights,
		     double *ans){
  /* This computes the least squares for the Brown-Resnick model */
  
  int i;
  double *rho, res;

  rho = (double *)R_alloc(*nPairs, sizeof(double));

  *ans = -brownResnick(dist, *nPairs, *range, *smooth, rho);
    
  if (*ans != 0.0)
    return;

  for (i=0;i<*nPairs;i++){
    res = 2 * pnorm(0.5 * rho[i], 0.0, 1.0, 1, 0) - extcoeff[i];
    *ans += res * res / (weights[i] * weights[i]);
  }

  return;
}
