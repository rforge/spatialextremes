#include "header.h"

void distance(double *coord, int *nDim, int *nSite,
	      double *dist){

  //This function computes the euclidean distance between each pair of
  //locations
  int i, j, k, currentPair = 0;

  for (i=0;i<(*nSite-1);i++){
    for (j=i+1;j<*nSite;j++){
      for (k=0;k<*nDim;k++)
	dist[currentPair] = dist[currentPair] + 
	  R_pow_di(coord[i + k * *nSite] -
		   coord[j + k * *nSite], 2);

      dist[currentPair] = sqrt(dist[currentPair]);
      currentPair++;
    }
  }
}

void distVecFct(double *coord, int *nSite, int *nDim, 
		double *distVec){

  //This function computes the distance vector between each pair of
  //locations
  const int nPair = *nSite * (*nSite - 1) / 2;
  int i, j, k, currentPair = -1;
  
  for (i=0;i<(*nSite-1);i++){
    for (j=i+1;j<*nSite;j++){
      currentPair++;
      for (k=0;k<*nDim;k++)
	distVec[k * nPair + currentPair] = coord[k * *nSite + j] -
	  coord[k * *nSite + i];

    }
  }
  return;
}
	
      

int gev2frech(double *data, int nObs, int nSite, double *locs, 
	      double *scales, double *shapes, double *jac,
	      double *frech){

  //This function transforms the GEV observations to unit Frechet ones
  //and computes the log of the jacobian of each transformation
  //When flag == 1, the GEV parameters are invalid.
  
  int i, j;
  const double eps = R_pow(DOUBLE_EPS, 0.3);
  
  for (i=0;i<nSite;i++){
    for (j=0;j<nObs;j++){
      frech[i * nObs + j] = (data[i * nObs + j] - locs[i])/ scales[i];
      
      if(fabs(shapes[i]) <= eps){
	jac[i * nObs + j] = frech[i * nObs + j] - log(scales[i]);
	frech[i * nObs + j] = exp(frech[i * nObs + j]);
      }
      
      else {
	frech[i * nObs + j] = 1 + shapes[i] * frech[i * nObs + j];
	
	if (frech[i * nObs + j] <= 0) {
	  return 1;
	}
	
	else{
	  jac[i * nObs + j] = (1/ shapes[i] -1) * 
	    log(frech[i * nObs + j]) - log(scales[i]);
	  frech[i * nObs + j] = R_pow(frech[i * nObs + j], 1/ shapes[i]);
	}
      }
    }
  }
  return 0;
}

int dsgnmat2Param(double *locdsgnmat, double *scaledsgnmat,
		  double *shapedsgnmat, double *loccoeff, 
		  double *scalecoeff, double *shapecoeff,
		  int nSite, int nloccoeff, int nscalecoeff,
		  int nshapecoeff, double *locs, double *scales,
		  double *shapes){

  int i, j;

  for (i=0;i<nSite;i++){
       
    locs[i] = 0.0;
    scales[i] = 0.0;
    shapes[i] = 0.0;
    
    for (j=0;j<nloccoeff;j++)
      locs[i] = locs[i] + loccoeff[j] * locdsgnmat[i + nSite * j];
    
    for (j=0;j<nscalecoeff;j++)
      scales[i] = scales[i] + scalecoeff[j] * scaledsgnmat[i + nSite * j];
    
    for (j=0;j<nshapecoeff;j++)
      shapes[i] = shapes[i] + shapecoeff[j] * shapedsgnmat[i + nSite * j];
    
    if ((scales[i]<=0) || (shapes[i] <= -1)){
      //printf("scales[%i] = %f\n", i, scales[i]);
      return 1;
    }
  }

  return 0;
}
  
void gev(double *prob, int *n, double *locs, double *scales, double *shapes,
	 double *quant){
  
  int i;
  
  for (i=0;i<*n;i++){
    
    if (shapes[i] == 0)
      quant[i] = locs[i] - scales[i] * log(-log(*prob));

    else
      quant[i] = locs[i] + scales[i] * R_pow(-log(*prob), shapes[i] - 1) /
	shapes[i];
  }
}

	
