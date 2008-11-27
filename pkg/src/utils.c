#include "header.h"

void distance(double *coord, int *nDim, int *nSite,
	      int *vec, double *dist){

  //This function computes either the euclidean distance or the
  //distance vector between each pair of locations
  const int nPair = *nSite * (*nSite - 1) / 2;
  int i, j, k, currentPair = 0;

  if (*vec){
    for (i=0;i<(*nSite-1);i++){
      for (j=i+1;j<*nSite;j++){
	for (k=0;k<*nDim;k++)
	  dist[k * nPair + currentPair] = coord[k * *nSite + j] -
	    coord[k * *nSite + i];
	
	currentPair++;
      }
    }
  }

  else{
    for (i=0;i<(*nSite-1);i++){
      for (j=i+1;j<*nSite;j++){
	for (k=0;k<*nDim;k++)
	  dist[currentPair] += R_pow_di(coord[i + k * *nSite] -
					coord[j + k * *nSite], 2);
	
	dist[currentPair] = sqrt(dist[currentPair]);
	currentPair++;
      }
    }
  }
} 

double gev2frech(double *data, int nObs, int nSite, double *locs, 
		 double *scales, double *shapes, double *jac,
		 double *frech){

  //This function transforms the GEV observations to unit Frechet ones
  //and computes the log of the jacobian of each transformation
  //When ans > 0.0, the GEV parameters are invalid.
  
  int i, j;
  double ans = 0.0;
  
  for (i=0;i<nSite;i++){
    for (j=0;j<nObs;j++){
      frech[i * nObs + j] = (data[i * nObs + j] - locs[i])/ scales[i];
      
      if(fabs(shapes[i]) <= 1e-6){
	shapes[i] = 0.0;
	jac[i * nObs + j] = frech[i * nObs + j] - log(scales[i]);
	frech[i * nObs + j] = exp(frech[i * nObs + j]);
      }
      
      else {
	frech[i * nObs + j] = 1 + shapes[i] * frech[i * nObs + j];
	
	if (frech[i * nObs + j] <= 0) {
	  //printf("1 + shape * (data - loc) <= 0!\n");
	  ans += R_pow_di(1.01 - frech[i * nObs + j], 2);
	  frech[i * nObs + j] = 1e-3;
	}
	
	jac[i * nObs + j] = (1/ shapes[i] -1) * 
	  log(frech[i * nObs + j]) - log(scales[i]);
	frech[i * nObs + j] = R_pow(frech[i * nObs + j], 1/ shapes[i]);
	
      }
    }
  }

  return ans;
}

double dsgnmat2Param(double *locdsgnmat, double *scaledsgnmat,
		     double *shapedsgnmat, double *loccoeff, 
		     double *scalecoeff, double *shapecoeff,
		     int nSite, int nloccoeff, int nscalecoeff,
		     int nshapecoeff, double *locs, double *scales,
		     double *shapes){

  //This function computes the GEV parameters from the design matrix
  //when ans > 0.0, the GEV parameters are invalid
  int i, j;
  double ans = 0.0;

  for (i=0;i<nSite;i++){
       
    locs[i] = 0.0;
    scales[i] = 0.0;
    shapes[i] = 0.0;
    
    for (j=0;j<nloccoeff;j++)
      locs[i] += loccoeff[j] * locdsgnmat[i + nSite * j];
    
    for (j=0;j<nscalecoeff;j++)
      scales[i] += scalecoeff[j] * scaledsgnmat[i + nSite * j];
    
    for (j=0;j<nshapecoeff;j++)
      shapes[i] += shapecoeff[j] * shapedsgnmat[i + nSite * j];
    
    if (scales[i]<=0){
      ans += R_pow_di(1.01 - scales[i], 2);
      scales[i] = 1.0;
    }

    if (shapes[i] <= -1){
      ans += R_pow_di(shapes[i] - 0.02, 2);
      shapes[i] = 0.0;
    }
  }

  return ans;
}
  
void gev(double *prob, int *n, double *locs, double *scales, double *shapes,
	 double *quant){

  //This function computes the GEV quantiles
  
  int i;
  
  for (i=0;i<*n;i++){
    
    if (scales[i] <= 0){
      quant[i] = R_NaReal;

    }

    if (shapes[i] == 0)
      quant[i] = locs[i] - scales[i] * log(-log(*prob));

    else
      quant[i] = locs[i] + scales[i] * (R_pow(-log(*prob), -shapes[i]) - 1) /
	shapes[i];
  }
}

	
void dsgnmat2Alpha(double *alphadsgnmat, double *alphacoeff, 
		   int nSite, int nalphacoeff, double *alphas){

  //This function computes the 'alpha' values from the design matrix
  //the 'alphas' are used in the schlatherind model
  //We use the expit function to ensure that the alphas always lie in
  //(0,1)
  int i, j;
  double ans = 0.0;

  for (i=0;i<nSite;i++){
       
    alphas[i] = 0.0;
        
    for (j=0;j<nalphacoeff;j++)
      alphas[i] += alphacoeff[j] * alphadsgnmat[i + nSite * j];

    //We use the expit function to ensure that the alphas lie in [0,1]
    alphas[i] = exp(alphas[i]) / (1 + exp(alphas[i]));
  
  }
    
  return;
}
