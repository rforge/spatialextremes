#include "header.h"

void schlatherfull(int *covmod, double *data, double *dist, int *nSite,
		   int *nObs, double *locs, double *scales, double *shapes,
		   double *scale, double *smooth, double *dns){
  //This is the schlater model with a Withle-Matern covariance
  //function. It computes the pairwise log-likelihood
  
  struct toFrech toFrechObj;
  struct covComp covCompObj;
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i;
  double *jac, *rho;
  
  jac = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));

  //Some preliminary steps: Valid points?
  for (i=0;i<*nSite;i++){
    if ((scales[i] <= 0) || (shapes[i]<= -1)){
      //printf("scale GEV param. are ill-defined!\n");
      *dns = -1.0e35;
      return;
    }
  }
  
  //Stage 0: Compute the covariance at each location
  switch (*covmod){
  case 1:
    covCompObj = whittleMatern(dist, nPairs, *scale, *smooth);
    break;
  case 2:
    covCompObj = cauchy(dist, nPairs, *scale, *smooth);
    break;
  case 3:
    covCompObj = powerExp(dist, nPairs, *scale, *smooth);
    break;
  }
  
  if (covCompObj.flag == 1){
    *dns = -1.0e35;
    return;
  }
  
  rho = covCompObj.vec;
    
  //Stage 1: Transformation to unit Frechet
  toFrechObj = gev2frech(data, *nObs, *nSite, locs, scales, shapes);
    
  if (toFrechObj.flag == 1){
    *dns = -1.0e35;
    return;
  }

  data = toFrechObj.frech;
  jac = toFrechObj.jac;
    
  //Stage 2: Bivariate density computations
  *dns = lplikschlather(data, rho, jac, *nObs, *nSite);

}

void schlatherlm(int *covmod, double *data, double *dist, int *nDim, int *nSite, int *nObs,
		 double *locdsgnmat, int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
		 double *shapedsgnmat, int *nshapecoeff, double *loccoeff, double *scalecoeff,
		 double *shapecoeff, double *scale, double *smooth,
		 double *dns){
  //This is the schlater model
  //The GEV parameters are defined using a polynomial response surface
  
  struct toFrech toFrechObj;
  struct covComp covCompObj;
  struct toParam gevParam;
  const int nPairs = *nSite * (*nSite - 1) / 2;
  double *jac, *rho, *locs, *scales, *shapes;
  //c1, c2 and c3 are usefull quantities
  
  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  
  //Stage 1: Compute the covariance at each location
  switch (*covmod){
  case 1:
    covCompObj = whittleMatern(dist, nPairs, *scale, *smooth);
    break;
  case 2:
    covCompObj = cauchy(dist, nPairs, *scale, *smooth);
    break;
  case 3:
    covCompObj = powerExp(dist, nPairs, *scale, *smooth);
    break;
  }
  
  if (covCompObj.flag == 1){
    *dns = -1.0e35;
    return;
  }

  rho = covCompObj.vec;

  //Stage 2: Compute the GEV parameters using the design matrix
  gevParam = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			   loccoeff, scalecoeff, shapecoeff, *nSite,
			   *nloccoeff, *nscalecoeff, *nshapecoeff);

  if (gevParam.flag == 1){
    *dns = -1.0e35;
    return;
  }

  locs = gevParam.locs;
  scales = gevParam.scales;
  shapes = gevParam.shapes;

  //Stage 3: Transformation to unit Frechet
  toFrechObj = gev2frech(data, *nObs, *nSite, locs, scales, shapes);
    
  if (toFrechObj.flag == 1){
    *dns = -1.0e35;
    return;
  }

  data = toFrechObj.frech;
  jac = toFrechObj.jac;
  
  //Stage 4: Bivariate density computations
  *dns = lplikschlather(data, rho, jac, *nObs, *nSite);
  
}

void schlatherdsgnmat(int *covmod, double *data, double *dist, int *nDim, int *nSite, int *nObs,
		      double *locdsgnmat, double *locpenmat, int *nloccoeff, int *npparloc,
		      double *locpenalty, double *scaledsgnmat, double *scalepenmat,
		      int *nscalecoeff, int *npparscale, double *scalepenalty, double *shapedsgnmat,
		      double *shapepenmat, int *nshapecoeff, int *npparshape, double *shapepenalty,
		      double *loccoeff, double *scalecoeff, double *shapecoeff, double *scale,
		      double *smooth, double *dns){
  //This is the schlater model
  //The GEV parameters are defined using a polynomial response surface
  
  struct toFrech toFrechObj;
  struct covComp covCompObj;
  struct toParam gevParam;
  const int nPairs = *nSite * (*nSite - 1) / 2;
  double *jac, *rho, *locs, *scales, *shapes;
  //c1, c2 and c3 are usefull quantities
  
  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  
  //Stage 1: Compute the covariance at each location
  switch (*covmod){
  case 1:
    covCompObj = whittleMatern(dist, nPairs, *scale, *smooth);
    break;
  case 2:
    covCompObj = cauchy(dist, nPairs, *scale, *smooth);
    break;
  case 3:
    covCompObj = powerExp(dist, nPairs, *scale, *smooth);
    break;
  }
  
  if (covCompObj.flag == 1){
    *dns = -1.0e35;
    return;
  }

  rho = covCompObj.vec;

  //Stage 2: Compute the GEV parameters using the design matrix
  gevParam = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			   loccoeff, scalecoeff, shapecoeff, *nSite,
			   *nloccoeff, *nscalecoeff, *nshapecoeff);

  if (gevParam.flag == 1){
    *dns = -1.0e35;
    return;
  }

  locs = gevParam.locs;
  scales = gevParam.scales;
  shapes = gevParam.shapes;

  //Stage 3: Transformation to unit Frechet
  toFrechObj = gev2frech(data, *nObs, *nSite, locs, scales, shapes);
    
  if (toFrechObj.flag == 1){
    *dns = -1.0e35;
    return;
  }

  data = toFrechObj.frech;
  jac = toFrechObj.jac;
  
  //Stage 4: Bivariate density computations
  *dns = lplikschlather(data, rho, jac, *nObs, *nSite);

  if (*dns == -1.0e35){
    //printf("problem with the pairwise lik.\n");
    return;
  }

  //Stage 5: Removing the penalizing terms (if any)
  // 1- For the location parameter
  if (*locpenalty > 0)
    *dns = *dns - penalization(locpenmat, loccoeff, *locpenalty,
			       *nloccoeff, *npparloc);
 
  // 2- For the scale parameter
  if (*scalepenalty > 0)    
    *dns = *dns - penalization(scalepenmat, scalecoeff, *scalepenalty,
			       *nscalecoeff, *npparscale);
    
  // 3- For the shape parameter
  if (*shapepenalty > 0)
    *dns = *dns - penalization(shapepenmat, shapecoeff, *shapepenalty,
			       *nshapecoeff, *npparshape);
  
}
