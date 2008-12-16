#include "header.h"

void schlatherindfull(int *covmod, double *data, double *dist, int *nSite,
		      int *nObs, double *locs, double *scales, double *shapes,
		      double *alpha, double *sill, double *range, double *smooth,
		      int *fitmarge,double *dns){
  //This is the independent Schlater's model. It's a wrapper to several
  //sub-functions. It's named xxxfull as it either assume that the
  //margins are unit Frechet, or the GEV parameters are estimated at
  //each locations.
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i;
  double *jac, *rho, *frech;
  
  jac = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));
  frech = (double *)R_alloc(*nSite * *nObs, sizeof(double));

  *dns = 1.0;
  //Some preliminary steps: Valid points?
  if (*fitmarge){
    for (i=0;i<*nSite;i++){
      if (scales[i] <= 0){
	//printf("scales <= 0!!!\n");
	*dns += R_pow_di(1 - scales[i], 2);
	scales[i] = 1;
      }
      
      if (shapes[i] <= -1){
	//printf("shapes <= -1!!!\n");
	*dns += R_pow_di(shapes[i], 2);
	shapes[i] = 0.0;
      }
    }
  }

  if ((*alpha < 0) || (*alpha > 1)){
    if (*alpha < 0){
      *dns += R_pow_di(1 - *alpha, 2);
      *alpha = 0.0;
    }

    else{
      *dns += R_pow_di(*alpha, 2);
      *alpha = 1.0;
    }
  }    
   
  //Stage 0: Compute the covariance at each location
  switch (*covmod){
  case 1:
    *dns += whittleMatern(dist, nPairs, *sill, *range, *smooth, rho);
    break;
  case 2:
    *dns += cauchy(dist, nPairs, *sill, *range, *smooth, rho);
    break;
  case 3:
    *dns += powerExp(dist, nPairs, *sill, *range, *smooth, rho);
    break;
  }
  
  //Stage 1: Transformation to unit Frechet
  if (*fitmarge){
    *dns += gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		      jac, frech);
    *dns *= lplikschlatherind(frech, *alpha, rho, jac, *nObs, *nSite);
  }
  
  else {
    for (i=0;i<(*nSite * *nObs);i++)
      jac[i] = 0.0;

    *dns *= lplikschlatherind(data, *alpha, rho, jac, *nObs, *nSite);
  }
  
  if (!R_FINITE(*dns))
    *dns = MINF;

  return;

}

void schlatherinddsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs,
			 double *locdsgnmat, double *locpenmat, int *nloccoeff, int *npparloc,
			 double *locpenalty, double *scaledsgnmat, double *scalepenmat,
			 int *nscalecoeff, int *npparscale, double *scalepenalty, double *shapedsgnmat,
			 double *shapepenmat, int *nshapecoeff, int *npparshape, double *shapepenalty,
			 double *loccoeff, double *scalecoeff, double *shapecoeff, double *alpha,
			 double *sill, double *range, double *smooth, double *dns){
  //This is the independent Schlater's model.
  //The GEV parameters are defined using a polynomial response surface
  //or p-splines.
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  double *jac, *rho, *locs, *scales, *shapes, *frech;
    
  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  
  *dns = 1.0;

  if (*alpha < 0){
    *dns += R_pow_di(1 - *alpha, 2);
    *alpha = 0.0;
  }

  if (*alpha > 1){
    *dns += R_pow_di(*alpha, 2);
    *alpha = 1.0;
  }

  //Stage 1: Compute the covariance at each location
  switch (*covmod){
  case 1:
    *dns += whittleMatern(dist, nPairs, *sill, *range, *smooth, rho);
    break;
  case 2:
    *dns += cauchy(dist, nPairs, *sill, *range, *smooth, rho);
    break;
  case 3:
    *dns += powerExp(dist, nPairs, *sill, *range, *smooth, rho);
    break;
  }
  
  //Stage 2: Compute the GEV parameters using the design matrix
  *dns += dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			loccoeff, scalecoeff, shapecoeff, *nSite,
			*nloccoeff, *nscalecoeff, *nshapecoeff,
			locs, scales, shapes);

  //Stage 3: Transformation to unit Frechet
  *dns -= gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		    jac, frech);

  *dns *= lplikschlatherind(frech, *alpha, rho, jac, *nObs, *nSite);
    
  //Stage 5: Removing the penalizing terms (if any)
  // 1- For the location parameter
  if (*locpenalty > 0)
    *dns -= penalization(locpenmat, loccoeff, *locpenalty,
			 *nloccoeff, *npparloc);
  
  // 2- For the scale parameter
  if (*scalepenalty > 0)    
    *dns -= penalization(scalepenmat, scalecoeff, *scalepenalty,
			 *nscalecoeff, *npparscale);
  
  // 3- For the shape parameter
  if (*shapepenalty > 0)
    *dns -= penalization(shapepenmat, shapecoeff, *shapepenalty,
			 *nshapecoeff, *npparshape);
  
  if (!R_FINITE(*dns))
    *dns = MINF;

  return;
  
}
