#include "header.h"

double whittleMatern(double *dist, int nPairs, double sill, double range,
		     double smooth, double *rho){

  //This function computes the whittle-matern covariance function
  //between each pair of locations.
  //When flag == 1, the whittle-matern parameters are ill-defined.
  
  int i;
  double ans = 0.0;

  //Some preliminary steps: Valid points?
  if (smooth <= EPS){
    //printf("dependence parameters are ill-defined!\n");
    ans = R_pow_di(1 - smooth, 2) * MINF;
  }

  if (range <= EPS){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1 - range, 2) * MINF;
  }

  if (sill <= EPS){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1 - sill, 2) * MINF;
  }
  
  if (sill > 1){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(sill, 2) * MINF;
  }
  
  if (smooth > 50){
    //Required because it could lead to infinite rho values
    //printf("smooth is too large!\n");
    ans += R_pow_di(smooth - 50, 2) * MINF;
  }

  if (ans != 0.0)
    return ans;
  
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
  //When flag == 1, the cauchy parameters are ill-defined.

  int i;
  double ans = 0.0;

  //Some preliminary steps: Valid points?
  if (smooth <= 0){
    //printf("dependence parameters are ill-defined!\n");
    ans = R_pow_di(1 - smooth, 2) * MINF;
  }

  if (range <= EPS){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1 - range, 2) * MINF;
  }

  if (sill <= EPS){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1 - sill, 2) * MINF;
  }
  
  if (sill > 1){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(sill, 2) * MINF;
  }

  if (ans != 0.0)
    return ans;
    
  for (i=0;i<nPairs;i++)
    rho[i] = sill * R_pow(1 + R_pow_di(dist[i] / range, 2), -smooth);
    
  return ans;
}

double powerExp(double *dist, int nPairs, double sill, double range,
		double smooth, double *rho){

  //This function computes the powered exponential covariance function
  //between each pair of locations.
  //When flag == 1, the powered exponential parameters are ill-defined.

  int i;
  double ans = 0.0;
  
  //Some preliminary steps: Valid points?
  if (smooth < 0){
    //printf("dependence parameters are ill-defined!\n");
    ans = R_pow_di(1 - smooth, 2) * MINF;
  }

  if (range <= EPS){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1 - range, 2) * MINF;
  }

  if (sill <= EPS){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(1 - sill, 2) * MINF;
  }
  
  if (sill > 1){
    //printf("dependence parameters are ill-defined!\n");
    ans += R_pow_di(sill, 2) * MINF;
  }
  
  if (smooth > 2){
    //Required because it could lead to infinite rho values
    //printf("smooth is too large!\n");
    ans += R_pow_di(smooth - 1, 2) * MINF;
  }

  if (ans != 0.0)
    return ans;
  
  for (i=0;i<nPairs;i++)
    rho[i] = sill * exp(-R_pow(dist[i] / range, smooth));
    
  return ans;
}

double mahalDistFct(double *distVec, int nPairs, double *cov11,
		    double *cov12, double *cov22, double *mahal){
  //This function computes the mahalanobis distance between each pair
  //of locations. Currently this function is only valid in 2D
  //When flag == 1, the covariance matrix and/or the mahalanobis
  //distance is ill-defined.
  
  int i;
  double det;

  det = *cov11 * *cov22 - R_pow_di(*cov12, 2);
  //We test if the covariance matrix is *not* nonnegative
  //definite e.g. all minor determinant are negative or 0
  if (det <= 0){
    //printf("Det of Sigma <= 0!\n");
    return R_pow_di(1 - det, 2) * MINF;
  }

  if (*cov11 <= 0){
    //printf("Cov11 is <=0!\n");
    return R_pow_di(1 - *cov11, 2) * MINF;
  }
  
  for (i=0;i<nPairs;i++){

    mahal[i] = (*cov11 *  R_pow_di(distVec[nPairs + i], 2) -
		2 * *cov12 * distVec[i] * distVec[nPairs + i] +
		*cov22 * R_pow_di(distVec[i], 2)) / det;
    
    mahal[i] = sqrt(mahal[i]);
    //We test if the Mahalanobis distance is singular.
    if (!R_FINITE(mahal[i])){
      //printf("sqrt(mahal) = %f\n", mahal[i]);
      return R_pow_di(1 - mahal[i], 2) * MINF;
    }
  }
  
  return 0.0;
}

double mahalDistFct3d(double *distVec, int nPairs, double *cov11,
		      double *cov12, double *cov13, double *cov22, 
		      double *cov23, double *cov33, double *mahal){
  //This function computes the mahalanobis distance between each pair
  //of locations. Currently this function is only valid in 3D
  //When flag == 1, the covariance matrix and/or the mahalanobis
  //distance is ill-defined.
  
  int i;
  double det, detMin;

  det = *cov11 * *cov22 * *cov33 - R_pow_di(*cov12, 2) * *cov33 -
    *cov11 * R_pow_di(*cov23, 2) + 2 * *cov12 * *cov13 * *cov23 -
    R_pow_di(*cov13, 2) * *cov22;
  detMin = *cov11 * *cov22 - R_pow_di(*cov12, 2);
  //We test if the covariance matrix is *not* nonnegative
  //definite e.g. all minor determinant are negative or 0
  if (det <= 0){
    //printf("Covariance matrice is singular!\n");
    return R_pow_di(1 - det, 2) * MINF;
  }

  if (*cov11 <= 0){
    //printf("Covariance matrice is singular!\n");
    return R_pow_di(1 - *cov11, 2) * MINF;
  }

  if (detMin <= 0){
    //printf("Covariance matrice is singular!\n");
    return R_pow_di(1 - detMin, 2) * MINF;
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

    mahal[i] = sqrt(mahal[i]);
    //We test if the Mahalanobis distance is singular.
    if (!R_FINITE(mahal[i]))
      return R_pow_di(1 - mahal[i], 2) * MINF;
  }
  
  return 0.0;
}
