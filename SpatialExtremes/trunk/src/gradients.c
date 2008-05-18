#include "header.h"

void smithgrad2(double *data, double *distVec, int *nSite,
		int *nObs, double *locdsgnmat, int *nloccoeff,
		double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		int *nshapecoeff, double *loccoeff, double *scalecoeff,
		double *shapecoeff, double *icov11, double *icov12,
		double *icov22, int *fitmarge, double *grad){

  //This is the Smith model. It computes the gradient of the pairwise log-likelihood
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i, j, k, currentPair = -1, flag;
  double c1, c2, dA, B, dB, C, dC, D, dD,
    *mahalDist, *locs, *scales, *shapes,
    jacCommon, *jac, *frech;
  //c1, c2 are useful quantities
  //A, B, C, D are part of the log-bivariate density
  //dB, dC and dD are their derivatives with respect to the
  //mahalanobis distance
  //jacCommon is the common part of all the jacobians - i.e. the common
  //part when deriving with respect to icov11, icov12 or icov22

  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));

  //Computing the Mahalanobis distance
  flag = mahalDistFct(distVec, nPairs, icov11, icov12,
		      icov22, mahalDist);

  //if (covCompObj.flag == 1){
  //  *dns = -1e50;
  //  return;
  //}

  //Compute the GEV parameters using the design matrix
  if (*fitmarge){
    
    flag = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			 loccoeff, scalecoeff, shapecoeff, *nSite,
			 *nloccoeff, *nscalecoeff, *nshapecoeff,
			 locs, scales, shapes);
    
    //if (flag == 1){
    //  *dns = -1e50;
    //  return;
    //}
    

    //Stage 1: Transformation to unit Frechet
    flag = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		     jac, frech);
    
    //if (flag == 1){
    //  *dns = -1e50;
    //  return;
    //}
    
  }
  
  //Stage 2: Gradient computations;
   for (k=0;k<*nObs;k++){
     currentPair = -1;
     for (i=0;i<(*nSite-1);i++){
       for (j=i+1;j<*nSite;j++){
	 
	 currentPair++;
	 
	 c1 = log(frech[k + j * *nObs] / frech[k + i * *nObs]) /
	   mahalDist[currentPair] + mahalDist[currentPair] / 2;
	 c2 = log(frech[k + i * *nObs] / frech[k + j * *nObs]) /
	   mahalDist[currentPair] + mahalDist[currentPair] / 2;
	 
	 //A = - pnorm(c1, 0., 1., 1, 0) / frech[k + i * *nObs] -
	 //  - pnorm(c2, 0., 1., 1, 0) / frech[k + j * *nObs];
	 B = - dnorm(c1, 0., 1., 0) / mahalDist[currentPair] /
	   frech[k + i * *nObs] / frech[k + j * *nObs] +
	   pnorm(c2, 0., 1., 1, 0) / R_pow_di(frech[k + j * *nObs], 2) +
	   dnorm(c2, 0., 1., 0) / mahalDist[currentPair] / 
	   R_pow_di(frech[k + j * *nObs], 2);
	 C = - dnorm(c2, 0., 1., 0) / mahalDist[currentPair] /
	   frech[k + i * *nObs] / frech[k + j * *nObs] +
	   pnorm(c1, 0., 1., 1, 0) / R_pow_di(frech[k + i * *nObs], 2) +
	   dnorm(c1, 0., 1., 0) / mahalDist[currentPair] / 
	   R_pow_di(frech[k + i * *nObs], 2);
	 D = c2 * dnorm(c1, 0., 1., 0) / frech[k + j * *nObs] /
	   R_pow_di(mahalDist[currentPair] * frech[k + i * *nObs], 2) +
	   c1 * dnorm(c2, 0., 1., 0) / frech[k + i * *nObs] /
	   R_pow_di(mahalDist[currentPair] * frech[k + j * *nObs], 2);
	 
	 dA = - c2 * dnorm(c1, 0., 1., 0) / frech[k + i * *nObs] /
	   mahalDist[currentPair] - c1 * dnorm(c2, 0., 1., 0) / 
	   frech[k + j * *nObs] / mahalDist[currentPair];
	 dB = (R_pow_di(c1, 2) - 1) * dnorm(c2, 0., 1., 0) /
	   R_pow_di(mahalDist[currentPair] * frech[k + j * *nObs], 2) +
	   (1 + c1 * c2 ) * dnorm(c1, 0., 1., 0) / frech[k + i * *nObs] /
	   frech[k + j * *nObs] / R_pow_di(mahalDist[currentPair], 2);
	 dC = (R_pow_di(c2, 2) - 1) * dnorm(c1, 0., 1., 0) /
	   R_pow_di(mahalDist[currentPair] * frech[k + i * *nObs], 2) +
	   (1 + c1 * c2) * dnorm(c2, 0., 1., 0) / frech[k + i * *nObs] /
	   frech[k + j * *nObs] / R_pow_di(mahalDist[currentPair], 2);
	 dD = (c1 - c1 * R_pow_di(c2, 2) - 2 * c2) * dnorm(c1, 0., 1., 0) / 
	   R_pow_di(mahalDist[currentPair], 3) /
	   R_pow_di(frech[k + i * *nObs], 2) / frech[k + j * *nObs] +
	   (c2 - R_pow_di(c1, 2) *c2 - 2 * c1) * dnorm(c2, 0., 1., 0) /
	   R_pow_di(mahalDist[currentPair], 3) / frech[k + i * *nObs] /
	   R_pow_di(frech[k + j * *nObs], 2);
	 
	 jacCommon = dA + (dB * C + B * dC + dD) / (B*C + D);
	 
	 grad[k] = grad[k] + R_pow_di(distVec[currentPair] , 2) * jacCommon;
	 grad[*nObs + k] = grad[*nObs + k] + 2 * distVec[currentPair] *
	   distVec[nPairs + currentPair] * jacCommon;
	 grad[2 * *nObs + k] = grad[2 * *nObs + k] + 
	   R_pow_di(distVec[nPairs + currentPair] , 2) * jacCommon;
       }
     }
   }
   
   return;
}
      
  
void smithgrad(double *data, double *distVec, int *nSite,
	       int *nObs, double *locdsgnmat, int *nloccoeff,
	       double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
	       int *nshapecoeff, double *loccoeff, double *scalecoeff,
	       double *shapecoeff, double *cov11, double *cov12,
	       double *cov22, int *fitmarge, double *grad){

  //This is the Smith model. It computes the gradient of the pairwise log-likelihood
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i, j, k, currentPair = -1, flag;
  double c1, c2, dA, B, dB, C, dC, D, dD,
    *mahalDist, *locs, *scales, *shapes,
    jacCommon, det, *jac, *frech;
  //c1, c2 are useful quantities
  //A, B, C, D are part of the log-bivariate density
  //dB, dC and dD are their derivatives with respect to the
  //mahalanobis distance
  //jacCommon is the common part of all the jacobians - i.e. the common
  //part when deriving with respect to cov11, cov12 or cov22
  //det is the determinant of the covariance matrix

  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));

  det = *cov11 * *cov22 - R_pow_di(*cov12, 2);

  //Computing the Mahalanobis distance
  flag = mahalDistFct(distVec, nPairs, cov11, cov12,
		      cov22, mahalDist);

  //if (flag == 1){
  //  *dns = -1e50;
  //  return;
  //}

  //Compute the GEV parameters using the design matrix
  if (*fitmarge){
    
    flag = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			 loccoeff, scalecoeff, shapecoeff, *nSite,
			 *nloccoeff, *nscalecoeff, *nshapecoeff,
			 locs, scales, shapes);
    
    //if (flag == 1){
    //  *dns = -1e50;
    //  return;
    //}
    
    //Stage 1: Transformation to unit Frechet
    flag = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		     jac, frech);
    
    //if (toFrechObj.flag == 1){
    //  *dns = -1e50;
    //  return;
    //}
    
  }
  
  //Stage 2: Gradient computations;
   for (k=0;k<*nObs;k++){
     currentPair = -1;
     for (i=0;i<(*nSite-1);i++){
       for (j=i+1;j<*nSite;j++){
	 
	 currentPair++;
	 
	 c1 = log(frech[k + j * *nObs] / frech[k + i * *nObs]) /
	   mahalDist[currentPair] + mahalDist[currentPair] / 2;
	 c2 = log(frech[k + i * *nObs] / frech[k + j * *nObs]) /
	   mahalDist[currentPair] + mahalDist[currentPair] / 2;
	 
	 //A = - pnorm(c1, 0., 1., 1, 0) / frech[k + i * *nObs] -
	 //  - pnorm(c2, 0., 1., 1, 0) / frech[k + j * *nObs];
	 B = - dnorm(c1, 0., 1., 0) / mahalDist[currentPair] /
	   frech[k + i * *nObs] / frech[k + j * *nObs] +
	   pnorm(c2, 0., 1., 1, 0) / R_pow_di(frech[k + j * *nObs], 2) +
	   dnorm(c2, 0., 1., 0) / mahalDist[currentPair] / 
	   R_pow_di(frech[k + j * *nObs], 2);
	 C = - dnorm(c2, 0., 1., 0) / mahalDist[currentPair] /
	   frech[k + i * *nObs] / frech[k + j * *nObs] +
	   pnorm(c1, 0., 1., 1, 0) / R_pow_di(frech[k + i * *nObs], 2) +
	   dnorm(c1, 0., 1., 0) / mahalDist[currentPair] / 
	   R_pow_di(frech[k + i * *nObs], 2);
	 D = c2 * dnorm(c1, 0., 1., 0) / frech[k + j * *nObs] /
	   R_pow_di(mahalDist[currentPair] * frech[k + i * *nObs], 2) +
	   c1 * dnorm(c2, 0., 1., 0) / frech[k + i * *nObs] /
	   R_pow_di(mahalDist[currentPair] * frech[k + j * *nObs], 2);
	 
	 dA = - c2 * dnorm(c1, 0., 1., 0) / frech[k + i * *nObs] /
	   mahalDist[currentPair] - c1 * dnorm(c2, 0., 1., 0) / 
	   frech[k + j * *nObs] / mahalDist[currentPair];
	 dB = (R_pow_di(c1, 2) - 1) * dnorm(c2, 0., 1., 0) /
	   R_pow_di(mahalDist[currentPair] * frech[k + j * *nObs], 2) +
	   (1 + c1 * c2 ) * dnorm(c1, 0., 1., 0) / frech[k + i * *nObs] /
	   frech[k + j * *nObs] / R_pow_di(mahalDist[currentPair], 2);
	 dC = (R_pow_di(c2, 2) - 1) * dnorm(c1, 0., 1., 0) /
	   R_pow_di(mahalDist[currentPair] * frech[k + i * *nObs], 2) +
	   (1 + c1 * c2) * dnorm(c2, 0., 1., 0) / frech[k + i * *nObs] /
	   frech[k + j * *nObs] / R_pow_di(mahalDist[currentPair], 2);
	 dD = (c1 - c1 * R_pow_di(c2, 2) - 2 * c2) * dnorm(c1, 0., 1., 0) / 
	   R_pow_di(mahalDist[currentPair], 3) /
	   R_pow_di(frech[k + i * *nObs], 2) / frech[k + j * *nObs] +
	   (c2 - R_pow_di(c1, 2) *c2 - 2 * c1) * dnorm(c2, 0., 1., 0) /
	   R_pow_di(mahalDist[currentPair], 3) / frech[k + i * *nObs] /
	   R_pow_di(frech[k + j * *nObs], 2);
	 
	 jacCommon = dA + (dB * C + B * dC + dD) / (B*C + D);
	 
	 grad[k] = grad[k] - R_pow_di((*cov11 * distVec[nPairs + currentPair] -
				       *cov22 * distVec[currentPair]) / det, 2) *
	   jacCommon;
	 grad[*nObs + k] = grad[*nObs + k] + 2 * 
	   (*cov11 * distVec[nPairs + currentPair] - *cov12 * distVec[currentPair]) *
	   (*cov12 * distVec[nPairs + currentPair] - *cov22 * distVec[currentPair]) /
	   R_pow_di(det, 2) * jacCommon;
	 grad[2 * *nObs + k] = grad[2 * *nObs + k] - 
	   R_pow_di((*cov11 * distVec[nPairs + currentPair] - *cov12 * distVec[currentPair]) /
		    det, 2) * jacCommon;
       }
     }
   }
   
   return;
}


void schlathergrad(int *covmod, double *data, double *dist, int *nSite,
		   int *nObs, double *locdsgnmat, int *nloccoeff,
		   double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		   int *nshapecoeff, double *loccoeff, double *scalecoeff,
		   double *shapecoeff, double *scale, double *smooth,
		   int *fitmarge, double *grad){

  //This is the Smith model. It computes the gradient of the pairwise log-likelihood
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i, j, k, currentPair = -1, flag;
  double c1, dA, B, dB, C, dC, D, dD,
    *rho, *locs, *scales, *shapes,
    jacCommon, *jac, *frech;
  //c1 is a useful quantity
  //A, B, C, D are part of the log-bivariate density
  //dB, dC and dD are their derivatives with respect to the
  //covariance function
  //jacCommon is the common part of all the jacobians - i.e. the common
  //part when deriving with respect to cov11, cov12 or cov22
  
  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));

  //Stage 0: Compute the covariance at each location
  switch (*covmod){
  case 1:
    flag = whittleMatern(dist, nPairs, *scale, *smooth, rho);
    break;
  case 2:
    flag = cauchy(dist, nPairs, *scale, *smooth, rho);
    break;
  case 3:
    flag = powerExp(dist, nPairs, *scale, *smooth, rho);
    break;
  }
  
  //Compute the GEV parameters using the design matrix
  if (*fitmarge){
    
    flag = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			 loccoeff, scalecoeff, shapecoeff, *nSite,
			 *nloccoeff, *nscalecoeff, *nshapecoeff,
			 locs, scales, shapes);
    
    //Stage 1: Transformation to unit Frechet
    flag = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		     jac, frech);
    
  }
  
  //Stage 2: Gradient computations;
   for (k=0;k<*nObs;k++){
     currentPair = -1;
     for (i=0;i<(*nSite-1);i++){
       for (j=i+1;j<*nSite;j++){
	 
	 currentPair++;
	 
	 c1 = sqrt(R_pow_di(frech[k + i * *nObs], 2) + 
		   R_pow_di(frech[k + j * *nObs], 2) -
		   2 * frech[k + i * *nObs] * frech[k + j * *nObs] *
		   rho[currentPair]);

	 B = (1 - R_pow_di(rho[currentPair], 2)) / 2 /
	   R_pow_di(c1, 3);
	 C = - (rho[currentPair] * frech[k + i * *nObs] - c1 -
		frech[k + j * *nObs]) / 2 / c1 /
	   R_pow_di(frech[k + i * *nObs], 2);
	 D = - (rho[currentPair] * frech[k + j * *nObs] - c1 -
		frech[k + i * *nObs]) / 2 / c1 /
	   R_pow_di(frech[k + j * *nObs], 2);

	 dA =  1 / 2 / c1;
	 dB = - rho[currentPair] / R_pow_di(c1, 3) + 3 * 
	   (1 - rho[currentPair]) * frech[k + i * *nObs] *
	   frech[k + j * *nObs] / R_pow_di(c1, 5);
	 dC = - (frech[k + i * *nObs] - frech[k + j * *nObs] *
		 rho[currentPair]) / 2 / R_pow_di(c1, 3);
	 dD = - (frech[k + j * *nObs] - frech[k + i * *nObs] *
		 rho[currentPair]) / 2 / R_pow_di(c1, 3);
		   

	 jacCommon = dA + (dB * C + B * dC + dD) / (B*C + D);
	 
	 switch (*covmod){
	 case 1:
	   grad[k] = grad[k] - rho[currentPair] * *smooth / *scale -
	     rho[currentPair] / R_pow_di(*scale, 2) + R_pow(2, 1 - *smooth) *
	     R_pow(dist[currentPair] / *scale, *smooth) * dist[currentPair] * 
	     bessel_k(dist[currentPair] / *scale, *smooth, 1) / 
	     R_pow_di(*scale, 2) * jacCommon;
	   
	   //I found no analytical form for the partial derivative of
	   //the Whittle-Matern covariance function with respect to
	   //the smooth parameter. Thus, I use a finite difference
	   //approach.
	   grad[*nObs + k] = grad[*nObs + k] - M_LN2 * rho[currentPair] -
	     rho[currentPair] * digamma(*smooth) + log(dist[currentPair] / *scale) *
	     rho[currentPair] + R_pow(2, 1 - *smooth) * R_pow(dist[currentPair] / *scale, *smooth) *
	     (bessel_k(dist[currentPair] / *scale, *smooth + 0.001, 1) - 
	      bessel_k(dist[currentPair] / *scale, *smooth, 1)) / 0.001 / gammafn(*smooth) * jacCommon;
	   break;
	 case 2:
	   grad[k] = grad[k] + rho[currentPair] * *smooth * 
	     R_pow_di(dist[currentPair] / *scale, 2) * R_pow(rho[currentPair], 1 / *smooth) *
	     jacCommon; 
	   grad[*nObs + k] = grad[*nObs + k] - rho[currentPair] * 
	     log(1 + R_pow_di(dist[currentPair], 2) / *scale) * jacCommon;
	   break;
	 case 3:
	   grad[k] = grad[k] + R_pow(dist[currentPair] / *scale, *smooth) *
	     *smooth * rho[currentPair] / *scale * jacCommon;
	   grad[*nObs + k] = grad[*nObs + k] - R_pow(dist[currentPair] / *scale, *smooth) *
	     log(dist[currentPair] / *scale) * rho[currentPair] * jacCommon;
	   break;	   
	 }
       }
     }
   }
   
   return;
}
      
  
