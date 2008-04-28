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
	
      

struct toFrech gev2frech(double *data, int nObs, int nSite, double *locs, 
			 double *scales, double *shapes){

  //This function transforms the GEV observations to unit Frechet ones
  //and computes the log of the jacobian of each transformation
  //When flag == 1, the GEV parameters are invalid.
  
  int i, j;
  struct toFrech ans;
  ans.frech =  (double *)R_alloc(nObs * nSite, sizeof(double));
  ans.jac =  (double *)R_alloc(nObs * nSite, sizeof(double));
  ans.flag = 0;  
  
  for (i=0;i<nSite;i++){
    for (j=0;j<nObs;j++){
      ans.frech[i * nObs + j] = (data[i * nObs + j] - locs[i])/ scales[i];
      
      if(shapes[i] == 0){
	ans.jac[i * nObs + j] = ans.frech[i * nObs + j] - log(scales[i]);
	ans.frech[i * nObs + j] = exp(ans.frech[i * nObs + j]);
      }
      
      else {
	ans.frech[i * nObs + j] = 1 + shapes[i] * ans.frech[i * nObs + j];
	
	if(ans.frech[i * nObs + j] <= 0) {
	  ans.flag = 1;
	  return ans;
	}
	
	else{
	  ans.jac[i * nObs + j] = (1/ shapes[i] -1) * 
	    log(ans.frech[i * nObs + j]) - log(scales[i]);
	  ans.frech[i * nObs + j] = R_pow(ans.frech[i * nObs + j], 1/ shapes[i]);
	}
      }
    }
  }
  return ans;
}

struct toParam dsgnmat2Param(double *locdsgnmat, double *scaledsgnmat,
			     double *shapedsgnmat, double *loccoeff, 
			     double *scalecoeff, double *shapecoeff,
			     int nSite, int nloccoeff, int nscalecoeff,
			     int nshapecoeff){

  int i, j;
  struct toParam ans;
  ans.locs =  (double *)R_alloc(nSite, sizeof(double));
  ans.scales =  (double *)R_alloc(nSite, sizeof(double));
  ans.shapes =  (double *)R_alloc(nSite, sizeof(double));
  ans.flag = 0;  
  
  for (i=0;i<nSite;i++){
    ans.locs[i] = 0.;
    ans.scales[i] = 0.;
    ans.shapes[i] = 0.;
    
    for (j=0;j<nloccoeff;j++)
      ans.locs[i] = ans.locs[i] + loccoeff[j] * locdsgnmat[i + nSite * j];
    
    for (j=0;j<nscalecoeff;j++)
      ans.scales[i] = ans.scales[i] + scalecoeff[j] * scaledsgnmat[i + nSite * j];
    
    for (j=0;j<nshapecoeff;j++)
      ans.shapes[i] = ans.shapes[i] + shapecoeff[j] * shapedsgnmat[i + nSite * j];
    
    if ((ans.scales[i]<=0) || (ans.shapes[i] <= -1)){
      //printf("scales[%i] = %f\n", i, scales[i]);
      ans.flag = 1;
      return ans;
    }
  }

  return(ans);
}
  
