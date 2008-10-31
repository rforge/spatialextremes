#include "header.h"

double whittleMatern(double *dist, int nPairs, double sill, double range,
		     double smooth, double *rho){

  //This function computes the whittle-matern covariance function
  //between each pair of locations.
  //When ans != 0.0, the whittle-matern parameters are ill-defined.
  
  int i;
  double ans = 0.0;

  //Some preliminary steps: Valid points?
  if (smooth < EPS){
    //printf("dependence parameters are ill-defined!\n");
    ans = R_pow_di(1 - smooth + EPS, 2);
    smooth = EPS;
  }

  if (range < EPS){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1 - range + EPS, 2);
    range = EPS;
  }

  if (sill < 1e-3){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1.001 - sill, 2);
    sill = 1e-3;
  }
  
  if (sill > 1){
    //printf("dependence parameters are ill-defined!\n");
    //the 1.02 factor is here to avoid problem with non
    //feasible region
    ans += R_pow_di(1.02 * sill, 2);
    sill = 1.0;
  }
  
  if (smooth > 150){
    //Required because it could lead to infinite rho values
    //printf("smooth is too large!\n");
    ans += R_pow_di(smooth - 150, 2);
    smooth = 150;
  }

  for (i=0;i<nPairs;i++){

    rho[i] = sill * R_pow(2, 1 - smooth) / gammafn(smooth) *
      R_pow(dist[i] / range, smooth) * 
      bessel_k(dist[i] / range, smooth, 1);
    
  }

  return ans;
}

double cauchy(double *dist, int nPairs, double sill, double range,
	      double smooth, double *rho){

  //This function computes the cauchy covariance function between each
  //pair of locations.
  //When ans != 0.0, the cauchy parameters are ill-defined.

  int i;
  double ans = 0.0;

  //Some preliminary steps: Valid points?
  if (smooth < 0){
    //printf("dependence parameters are ill-defined!\n");
    ans = R_pow_di(1 - smooth, 2);
    smooth = 1e-3;
  }

  if (range < EPS){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1 - range + EPS, 2);
    range = EPS;
  }

  if (sill < 1e-3){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1.001 - sill, 2);
    sill = 1e-3;
  }
  
  if (sill > 1){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1.02 * sill, 2);
    sill = 1.0;
  }

  for (i=0;i<nPairs;i++)
    rho[i] = sill * R_pow(1 + R_pow_di(dist[i] / range, 2), -smooth);
    
  return ans;
}

double powerExp(double *dist, int nPairs, double sill, double range,
		double smooth, double *rho){

  //This function computes the powered exponential covariance function
  //between each pair of locations.
  //When ans != 0.0, the powered exponential parameters are ill-defined.

  int i;
  double ans = 0.0;
  
  //Some preliminary steps: Valid points?
  if (smooth < 0){
    //printf("dependence parameters are ill-defined!\n");
    ans = R_pow_di(1 - smooth, 2);
    smooth = 1e-3;
  }

  if (range < EPS){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1 - range + EPS, 2);
    range = EPS;
  }

  if (sill < 1e-3){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1.001 - sill, 2);
    sill = 1e-3;
  }
  
  if (sill > 1){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1.02 * sill, 2);
    sill = 1.0;
  }
  
  if (smooth > 2){
    //Required because it could lead to infinite rho values
    //printf("smooth is too large!\n");
    ans += R_pow_di(smooth - 1, 2);
    smooth = 2.0;
  }

  for (i=0;i<nPairs;i++)
    rho[i] = sill * exp(-R_pow(dist[i] / range, smooth));
    
  return ans;
}

double mahalDistFct(double *distVec, int nPairs, double *cov11,
		    double *cov12, double *cov22, double *mahal){
  //This function computes the mahalanobis distance between each pair
  //of locations. Currently this function is only valid in 2D
  //When ans != 0.0, the covariance matrix and/or the mahalanobis
  //distance is ill-defined.
  
  int i;
  double det, ans = 0.0;

  det = *cov11 * *cov22 - R_pow_di(*cov12, 2);
  //We test if the covariance matrix is *not* nonnegative
  //definite e.g. all minor determinant are negative or 0
  if (det <= 1e-10){
    //printf("Det of Sigma <= 0!\n");
    ans += R_pow_di(1 - det + 1e-10, 2);
  }

  if (*cov11 <= 0){
    //printf("Cov11 is <=0!\n");
    ans += R_pow_di(1 - *cov11, 2);
  }
  
  for (i=0;i<nPairs;i++){

    mahal[i] = (*cov11 *  R_pow_di(distVec[nPairs + i], 2) -
		2 * *cov12 * distVec[i] * distVec[nPairs + i] +
		*cov22 * R_pow_di(distVec[i], 2)) / det;
    
    //We test if the Mahalanobis distance is singular.
    if (!R_FINITE(sqrt(mahal[i]))){
      //printf("sqrt(mahal) = %f\n", mahal[i]);
      ans += R_pow_di(1 + fabs(mahal[i]), 2);
      mahal[i] = EPS;
    }

    mahal[i] = sqrt(mahal[i]);
  }
  
  return ans;
}

double mahalDistFct3d(double *distVec, int nPairs, double *cov11,
		      double *cov12, double *cov13, double *cov22, 
		      double *cov23, double *cov33, double *mahal){
  //This function computes the mahalanobis distance between each pair
  //of locations. Currently this function is only valid in 3D
  //When ans != 0.0, the covariance matrix and/or the mahalanobis
  //distance is ill-defined.
  
  int i;
  double det, detMin, ans = 0.0;

  det = *cov11 * *cov22 * *cov33 - R_pow_di(*cov12, 2) * *cov33 -
    *cov11 * R_pow_di(*cov23, 2) + 2 * *cov12 * *cov13 * *cov23 -
    R_pow_di(*cov13, 2) * *cov22;
  detMin = *cov11 * *cov22 - R_pow_di(*cov12, 2);
  //We test if the covariance matrix is *not* nonnegative
  //definite e.g. all minor determinant are negative or 0
  if (det <= 1e-10){
    //printf("Covariance matrice is singular!\n");
    ans += R_pow_di(1 - det + 1e-10, 2);
  }

  if (*cov11 <= 0){
    //printf("Covariance matrice is singular!\n");
    ans += R_pow_di(1 - *cov11, 2);
  }

  if (detMin <= 0){
    //printf("Covariance matrice is singular!\n");
    ans += R_pow_di(1 - detMin, 2);
  }
  
  for (i=0;i<nPairs;i++){

    mahal[i] = (*cov11 * *cov22 * R_pow_di(distVec[2 * nPairs + i], 2) -
		R_pow_di(*cov12 * distVec[2 * nPairs + i], 2) - 2 * *cov11 *
		*cov23 * distVec[nPairs + i] * distVec[2 * nPairs + i] + 2 *
		*cov12 * *cov13 * distVec[nPairs + i] * distVec[2 * nPairs + i] +
		2 * *cov12 * *cov23 * distVec[i] * distVec[2 * nPairs + i] - 2 *
		*cov13 * *cov22 * distVec[i] * distVec[2 * nPairs + i] + *cov11 *
		*cov33 * R_pow_di(distVec[nPairs + i], 2) - 
		R_pow_di(*cov13 * distVec[nPairs + i], 2) - 2 * *cov12 * *cov33 *
		distVec[i] * distVec[nPairs + i] + 2 * *cov13 * *cov23 *
		distVec[i] * distVec[nPairs + i] + *cov22 * *cov33 * 
		R_pow_di(distVec[i], 2) - R_pow_di(*cov23 * distVec[i], 2)) / det;

    //We test if the Mahalanobis distance is singular.
    if (!R_FINITE(sqrt(mahal[i]))){
      ans += R_pow_di(1 + fabs(mahal[i]), 2);
      mahal[i] = EPS;
    }

    mahal[i] = sqrt(mahal[i]);
  }
  
  return ans;
}

double geomCovariance(double *dist, int nPairs, int covmod,
		      double sigma2, double sill, double range,
		      double smooth, double *rho){
  int i;
  double ans;

  switch (covmod){
  case 1:
    ans = whittleMatern(dist, nPairs, sill, range, smooth, rho);
    break;
  case 2:
    ans = cauchy(dist, nPairs, sill, range, smooth, rho);
    break;
  case 3:
    ans = powerExp(dist, nPairs, sill, range, smooth, rho);
    break;
  }

  if (sigma2 > 3){
    ans += R_pow_di(ans-2, 2);
    ans = 3.0;
  }

  else
    for (i=0;i<nPairs;i++)
      rho[i] = 2 * sigma2 * (1 - rho[i]);

  return ans;
}

      

  
  
