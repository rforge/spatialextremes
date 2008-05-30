#include "header.h"

void smithfull(double *data, double *distVec, int *nSite,
	       int *nObs, double *locs, double *scales, double *shapes,
	       double *icov11, double *icov12, double *icov22, double *dns){
  //This is the Smith model. It computes the pairwise log-likelihood

  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i, flag = 0;
  double *jac, *mahalDist, *frech;
  
  jac = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  frech = (double *)R_alloc(*nSite * *nObs, sizeof(double));

  //Some preliminary steps: Valid points?
  for (i=0;i<*nSite;i++)
    if ((scales[i] <= 0) || (shapes[i] <= -1)){
      //printf("scales <= 0!!!\n");
      *dns = -1.0e35;
      return;
    }

  //Stage 1: Computing the Mahalanobis distance
  flag = mahalDistFct(distVec, nPairs, icov11,
		      icov12, icov22, mahalDist);
  
  if (flag == 1){
    //printf("Problem with mahal. dist\n");
    *dns = -1.0e35;
    return;
  }

  //Stage 2: Transformation to unit Frechet
  flag = gev2frech(data, *nObs, *nSite, locs, scales,
		   shapes, jac, frech);
    
  if (flag == 1){
    //printf("problem with conversion to unit frechet\n");
    *dns = -1.0e35;
    return;
  }
  
  //Stage 3: Bivariate density computations
  *dns = lpliksmith(frech, mahalDist, jac, *nObs, *nSite);

  if (!R_FINITE(*dns))
    *dns = -1.0e35;

  return;
}

void smithdsgnmat(double *data, double *distVec, int *nSite, int *nObs, 
		  double *locdsgnmat, double *locpenmat, int *nloccoeff,
		  int *npparloc, double *locpenalty, double *scaledsgnmat,
		  double *scalepenmat, int *nscalecoeff, int *npparscale,
		  double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		  int *nshapecoeff, int *npparshape, double *shapepenalty,
		  double *loccoeff, double *scalecoeff, double *shapecoeff,
		  double *icov11, double *icov12, double *icov22, double *dns){
  //This is the Smith model. It computes the pairwise log-likelihood
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i, j, flag = 0;
  double *jac, *mahalDist, *locs, *scales, *shapes, *frech;
  
  jac = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  
  //Stage 1: Computing the Mahalanobis distance
  flag = mahalDistFct(distVec, nPairs, icov11, icov12,
		      icov22, mahalDist);

  if (flag == 1){
    //printf("problem with mahal. dist\n");
    *dns = -1.0e35;
    return;
  }

  //Stage 2: Computing the GEV parameters using the design matrix
  flag = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
		       loccoeff, scalecoeff, shapecoeff, *nSite,
		       *nloccoeff, *nscalecoeff, *nshapecoeff,
		       locs, scales, shapes);

  if (flag == 1){
    //printf("problem with the GEV parameters\n");
    *dns = -1.0e35;
    return;
  }

  //Stage 3: Transformation to unit Frechet
  flag = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		   jac, frech);
    
  if (flag == 1){
    //printf("problem with conversion to unit frechet\n");
    *dns = -1.0e35;
    return;
  }
  
  //Stage 4: Bivariate density computations
  *dns = lpliksmith(frech, mahalDist, jac, *nObs, *nSite);
  
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

  if (!R_FINITE(*dns))
    *dns = -1.0e35;

  return;
}
