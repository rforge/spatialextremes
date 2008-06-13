#include "header.h"

int whittleMatern(double *dist, int nPairs, double scale,
		   double smooth, double *rho){

  //This function computes the whittle-matern covariance function
  //between each pair of locations.
  //When flag == 1, the whittle-matern parameters are ill-defined.
  
  int i;

  //Some preliminary steps: Valid points?
  if ((smooth <= 0) || (scale <= 0)){
    //printf("dependence parameters are ill-defined!\n");
    return 1;
  }
  
  if (smooth > 171){
    //Required because gammafn is infinite
    //printf("smooth is too large!\n");
    return 1;
  }
  
  for (i=0;i<nPairs;i++){

    rho[i] = R_pow(2, 1 - smooth) / gammafn(smooth) *
      R_pow(dist[i] / scale, smooth) * 
      bessel_k(dist[i] / scale, smooth, 1);


    if (!R_FINITE(rho[i]))
      rho[i] = 1.0;

  }

  return 0;
}

int cauchy(double *dist, int nPairs, double scale,
	   double smooth, double *rho){

  //This function computes the cauchy covariance function between each
  //pair of locations.
  //When flag == 1, the cauchy parameters are ill-defined.

  int i;

  //Some preliminary steps: Valid points?
  if ((smooth <= 0) || (scale <= 0)){
    //printf("dependence parameters are ill-defined!\n");
    return 1;
  }
  
  for (i=0;i<nPairs;i++)
    rho[i] = R_pow(1 + R_pow_di(dist[i] / scale, 2), -smooth);
    
  return 0;
}

int powerExp(double *dist, int nPairs, double scale,
	     double smooth, double *rho){

  //This function computes the powered exponential covariance function
  //between each pair of locations.
  //When flag == 1, the powered exponential parameters are ill-defined.

  int i;
  
  //Some preliminary steps: Valid points?
  if ((smooth < 0) || (smooth > 2) || (scale <= 0)){
    //printf("dependence parameters are ill-defined!\n");
    return 1;
  }
  
  for (i=0;i<nPairs;i++)
    rho[i] = exp(-R_pow(dist[i] / scale, smooth));
    
  return 0;
}

int genHyper(double *dist, int nPairs, double scale,
	     double smooth1, double smooth2,
	     double smooth3, double *rho){

  //This function computes the generalized hyperbolic covariance function
  //between each pair of locations.
  //When flag == 1, the generalized hyperbolic parameters are ill-defined.

  int i;

  //Some preliminary steps: Valid points?
  if (((smooth3 > 0) && ((smooth2 <= 0) || (smooth1 < 0))) ||
      ((smooth3 == 0) && ((smooth2 <= 0) || (smooth1 <= 0))) ||
      ((smooth3 < 0) && ((smooth2 < 0) || (smooth1 <= 0)))){
    //printf("dependence parameters are ill-defined!\n");
    return 1;
  }
  
  for (i=0;i<nPairs;i++)
    rho[i] = sqrt(smooth2 / 2 / M_PI) / R_pow(smooth1, smooth3) /
      bessel_k(smooth1 * smooth2, smooth3, 1) * 
      R_pow(pythag(smooth1, dist[i]), smooth3 - .5) *
      bessel_k(smooth2 * pythag(smooth1, dist[i]), smooth3 - .5, 1);

  return 0;
}

int mahalDistFct(double *distVec, int nPairs, double *cov11,
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
  if ((det <= 0) || (*cov11 <= 0)){
    //printf("Covariance matrice is singular!\n");
    return 1;
  }
  
  for (i=0;i<nPairs;i++){

    mahal[i] = (*cov11 *  R_pow_di(distVec[nPairs + i], 2) -
		2 * *cov12 * distVec[i] * distVec[nPairs + i] +
		*cov22 * R_pow_di(distVec[i], 2)) / det;
    
    //We test if the Mahalanobis distance is singular.
    if (mahal[i] <= 0)
      return 1;
    
    mahal[i] = sqrt(mahal[i]);
  }
  
  return 0;
}

int mahalDistFct3d(double *distVec, int nPairs, double *cov11,
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
  if ((det <= 0) || (*cov11 <= 0) || (detMin <= 0)){
    //printf("Covariance matrice is singular!\n");
    return 1;
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
    if (mahal[i] <= 0)
      return 1;
    
    mahal[i] = sqrt(mahal[i]);
  }
  
  return 0;
}
