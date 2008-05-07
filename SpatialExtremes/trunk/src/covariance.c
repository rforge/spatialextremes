#include "header.h"

struct covComp whittleMatern(double *dist, int nPairs, double scale,
			     double smooth){

  //This function computes the whittle-matern covariance function
  //between each pair of locations.
  //When flag == 1, the whittle-matern parameters are ill-defined.
  
  int i;
  struct covComp ans;  
  ans.vec =  (double *)R_alloc(nPairs, sizeof(double));
  ans.flag = 0;  

  //Some preliminary steps: Valid points?
  if ((smooth <= 0) || (scale <= 0)){
    //printf("dependence parameters are ill-defined!\n");
    ans.flag = 1;
    return ans;
  }
  
  if (smooth > 171){
    //Required because gammafn is infinite
    //printf("smooth is too large!\n");
    ans.flag = 1;
    return ans;
  }
  
  for (i=0;i<nPairs;i++){

    ans.vec[i] = R_pow(2, 1 - smooth) / gammafn(smooth) *
	  R_pow(dist[i] / scale, smooth) * 
	  bessel_k(dist[i] / scale, smooth, 1);


    if (!R_FINITE(ans.vec[i]))
      ans.vec[i] = 1.0;

  }

  return ans;
}

struct covComp cauchy(double *dist, int nPairs, double scale,
		      double smooth){

  //This function computes the cauchy covariance function between each
  //pair of locations.
  //When flag == 1, the cauchy parameters are ill-defined.

  int i;
  struct covComp ans;  
  ans.vec =  (double *)R_alloc(nPairs, sizeof(double));
  ans.flag = 0;  

  //Some preliminary steps: Valid points?
  if ((smooth <= 0) || (scale <= 0)){
    //printf("dependence parameters are ill-defined!\n");
    ans.flag = 1;
    return ans;
  }
  
  for (i=0;i<nPairs;i++)
    ans.vec[i] = R_pow(1 + R_pow_di(dist[i] / scale, 2), -smooth);
    
  return ans;
}

struct covComp powerExp(double *dist, int nPairs, double scale,
			double smooth){

  //This function computes the powered exponential covariance function
  //between each pair of locations.
  //When flag == 1, the powered exponential parameters are ill-defined.

  int i;
  struct covComp ans;  
  ans.vec =  (double *)R_alloc(nPairs, sizeof(double));
  ans.flag = 0;  

  //Some preliminary steps: Valid points?
  if ((smooth < 0) || (smooth > 2) || (scale <= 0)){
    //printf("dependence parameters are ill-defined!\n");
    ans.flag = 1;
    return ans;
  }
  
  for (i=0;i<nPairs;i++)
    ans.vec[i] = exp(-R_pow(dist[i] / scale, smooth));
    
  return ans;
}

struct covComp genHyper(double *dist, int nPairs, double scale,
			double smooth1, double smooth2,
			double smooth3){

  //This function computes the generalized hyperbolic covariance function
  //between each pair of locations.
  //When flag == 1, the generalized hyperbolic parameters are ill-defined.

  int i;
  struct covComp ans;  
  ans.vec =  (double *)R_alloc(nPairs, sizeof(double));
  ans.flag = 0;  

  //Some preliminary steps: Valid points?
  if (((smooth3 > 0) && ((smooth2 <= 0) || (smooth1 < 0))) ||
      ((smooth3 == 0) && ((smooth2 <= 0) || (smooth1 <= 0))) ||
      ((smooth3 < 0) && ((smooth2 < 0) || (smooth1 <= 0)))){
    //printf("dependence parameters are ill-defined!\n");
    ans.flag = 1;
    return ans;
  }
  
  for (i=0;i<nPairs;i++)
    ans.vec[i] = sqrt(smooth2 / 2 / M_PI) / R_pow(smooth1, smooth3) /
      bessel_k(smooth1 * smooth2, smooth3, 1) * 
      R_pow(pythag(smooth1, dist[i]), smooth3 - .5) *
      bessel_k(smooth2 * pythag(smooth1, dist[i]), smooth3 - .5, 1);

  return ans;
}


struct covComp mahalDistFct2(double *distVec, int nPairs, double *icov11,
			    double *icov12, double *icov22){
  
  //This function computes the mahalanobis distance between each pair
  //of locations. Currently this function is only valid in 2D
  //When flag == 1, the covariance matrix and/or the mahalanobis
  //distance is ill-defined.
  int i;
  double det;
  struct covComp ans;

  ans.vec =  (double *)R_alloc(nPairs, sizeof(double));
  ans.flag = 0;  

  det = *icov11 * *icov22 - R_pow_di(*icov12, 2);
  //We test if the (inverse) covariance matrix is *not* nonnegative
  //definite e.g. all minor determinants are negative or 0
  if ((det <= 0) || (*icov22 <= 0)){
    ans.flag = 1;
    //printf("Covariance matrice is singular!\n");
    return ans;
  }
  
  for (i=0;i<nPairs;i++){

    ans.vec[i] = *icov11 * R_pow_di(distVec[i], 2) +
      2 * *icov12 * distVec[i] * distVec[nPairs + i] +
      *icov22 * R_pow_di(distVec[nPairs + i], 2);
    //We test if the Mahalanobis distance is singular.
    if (!R_FINITE(ans.vec[i]) || (ans.vec[i] <= 0)){
      ans.flag = 1;
      //printf("mahalDist^2 is erradic!\n"); 
      return ans;
    }
    
    ans.vec[i] = sqrt(ans.vec[i]);
  }
  
  return ans;
}

struct covComp mahalDistFct(double *distVec, int nPairs, double *cov11,
			     double *cov12, double *cov22){
  //This function computes the mahalanobis distance between each pair
  //of locations. Currently this function is only valid in 2D
  //When flag == 1, the covariance matrix and/or the mahalanobis
  //distance is ill-defined.
  //The difference from the previous function is that it uses the
  //covariance matrix instead of the inverse covariance matrix directly

  int i;
  const double eps = R_pow(DOUBLE_EPS, 0.3);
  double det;
  struct covComp ans;

  ans.vec =  (double *)R_alloc(nPairs, sizeof(double));
  ans.flag = 0;  

  det = *cov11 * *cov22 - R_pow_di(*cov12, 2);
  //We test if the covariance matrix is *not* nonnegative
  //definite e.g. all minor determinant are negative or 0
  if ((det <= eps) || (*cov22 <= 0)){
    ans.flag = 1;
    //printf("Covariance matrice is singular!\n");
    return ans;
  }
  
  for (i=0;i<nPairs;i++){

    ans.vec[i] = (*cov11 *  R_pow_di(distVec[nPairs + i], 2) -
			    2 * *cov12 * distVec[i] * distVec[nPairs + i] +
			    *cov22 * R_pow_di(distVec[i], 2)) / det;
    
    //We test if the Mahalanobis distance is singular.
    //if (!R_FINITE(ans.vec[i]) || (ans.vec[i] <= 0)){
    //  ans.flag = 1;
    //  //printf("mahalDist^2 is erradic!\n"); 
    //  return ans;
    //}
    
    ans.vec[i] = sqrt(ans.vec[i]);
  }
  
  return ans;
}
