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
	  dist[currentPair] += (coord[i + k * *nSite] -	coord[j + k * *nSite]) *
	    (coord[i + k * *nSite] - coord[j + k * *nSite]);
	
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
    
  for (i=nSite;i--;){
    if (shapes[i] == 0.0){
            
      for (j=nObs;j--;){
	frech[i * nObs + j] = (data[i * nObs + j] - locs[i])/ scales[i];
	jac[i * nObs + j] = frech[i * nObs + j] - log(scales[i]);
	frech[i * nObs + j] = exp(frech[i * nObs + j]);
      }
    }
      
    else {
      for (j=nObs;j--;){
	frech[i * nObs + j] = 1 + shapes[i] * (data[i * nObs + j] - locs[i]) /
	  scales[i];
	
	if (frech[i * nObs + j] <= 0)
	  return MINF;
	
	jac[i * nObs + j] = (1/ shapes[i] - 1) * 
	  log(frech[i * nObs + j]) - log(scales[i]);
	frech[i * nObs + j] = R_pow(frech[i * nObs + j], 1/ shapes[i]);
	
      }
    }
  }
  
  return 0.0;
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
    
    if ((scales[i]<=0) || (shapes[i] <= -1))
      return MINF;
  }

  return 0.0;
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

  for (i=0;i<nSite;i++){
       
    alphas[i] = 0.0;
        
    for (j=0;j<nalphacoeff;j++)
      alphas[i] += alphacoeff[j] * alphadsgnmat[i + nSite * j];

    //We use the expit function to ensure that the alphas lie in [0,1]
    alphas[i] = exp(alphas[i]) / (1 + exp(alphas[i]));
  
  }
    
  return;
}

void dsgnmat2Sigma2(double *sigma2dsgnmat, double *sigma2coeff, 
		    int nSite, int nsigma2coeff, double *sigma2){

  //This function computes the 'sigma2' values from the design matrix
  //the 'sigma2' are used in the non-stationary geometric gaussian model
  int i, j;
  
  for (i=0;i<nSite;i++){
       
    sigma2[i] = 0.0;
        
    for (j=0;j<nsigma2coeff;j++)
      sigma2[i] += sigma2coeff[j] * sigma2dsgnmat[i + nSite * j];

    //We use a log link function to ensure that the sigma2s lie are positive
    sigma2[i] = exp(sigma2[i]);
  
  }
    
  return;
}
