#include "header.h"

void smithfull(double *data, double *distVec, int *nSite, int *nObs,
	       int *weighted, double *weights, double *locs, double *scales,
	       double *shapes, double *cov11, double *cov12, double *cov22,
	       int *fitmarge, double *dns){
  //This is the Smith's model. It's a wrapper to several
  //sub-functions. It's named xxxfull as it either assume that the
  //margins are unit Frechet, or the GEV parameters are estimated at
  //each locations.

  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i;
  double *jac, *mahalDist, *frech;
  
  jac = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  frech = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  
  //Some preliminary steps: Valid points?
  if (*fitmarge){
    for (i=0;i<*nSite;i++){
      if ((scales[i] <= 0) || (shapes[i] <= -1)){
	*dns = MINF;
	return;
      }
    }
  }
  
  //Stage 1: Computing the Mahalanobis distance
  *dns = mahalDistFct(distVec, nPairs, cov11, cov12, cov22, mahalDist);

  if (*dns != 0.0)
      return;

  //Stage 2: Transformation to unit Frechet
  if (*fitmarge){
    *dns = gev2frech(data, *nObs, *nSite, locs, scales, shapes, jac, frech);

    if (*dns != 0.0)
      return;

    if (*weighted)
      *dns = wlpliksmith(frech, mahalDist, jac, *nObs, *nSite, weights);

    else
      *dns = lpliksmith(frech, mahalDist, jac, *nObs, *nSite);
  }
  
  else {
    memset(jac, 0, *nSite * *nObs * sizeof(double));
    
    if (*weighted)
      *dns = wlpliksmith(data, mahalDist, jac, *nObs, *nSite, weights);

    else
      *dns = lpliksmith(data, mahalDist, jac, *nObs, *nSite);
  }
  
  return;
}

void smithdsgnmat(double *data, double *distVec, int *nSite, int *nObs, int *weighted,
		  double *weights, double *locdsgnmat, double *locpenmat, int *nloccoeff,
		  int *npparloc, double *locpenalty, double *scaledsgnmat,
		  double *scalepenmat, int *nscalecoeff, int *npparscale,
		  double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		  int *nshapecoeff, int *npparshape, double *shapepenalty,
		  int *usetempcov, double *tempdsgnmatloc, double *temppenmatloc,
		  int *ntempcoeffloc, int *nppartempcoeffloc, double *temppenaltyloc,
		  double *tempdsgnmatscale, double *temppenmatscale, int *ntempcoeffscale,
		  int *nppartempcoeffscale, double *temppenaltyscale, double *tempdsgnmatshape,
		  double *temppenmatshape, int *ntempcoeffshape, int *nppartempcoeffshape,
		  double *temppenaltyshape, double *loccoeff, double *scalecoeff,
		  double *shapecoeff, double *tempcoeffloc, double *tempcoeffscale,
		  double *tempcoeffshape, double *cov11, double *cov12, double *cov22,
		  double *dns){
  //This is the Smith's model. It's named xxxdsgnmat as either linear
  //models or p-splines are used for the gev parameters.
  
  const int nPairs = *nSite * (*nSite - 1) / 2,
    flag = usetempcov[0] + usetempcov[1] + usetempcov[3];
  double *jac, *mahalDist, *locs, *scales, *shapes, *frech, *trendlocs, *trendscales,
    *trendshapes;
  
  jac = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  
  //Stage 1: Computing the Mahalanobis distance
  *dns = mahalDistFct(distVec, nPairs, cov11, cov12, cov22, mahalDist);

  if (*dns != 0)
    return;

  //Stage 2: Computing the GEV parameters using the design matrix
  *dns = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat, loccoeff, scalecoeff, shapecoeff,
		       *nSite, *nloccoeff, *nscalecoeff, *nshapecoeff, locs, scales, shapes);

  if (flag){
    int i, j;
    trendlocs = (double *)R_alloc(*nObs, sizeof(double));
    trendscales = (double *)R_alloc(*nObs, sizeof(double));
    trendshapes = (double *)R_alloc(*nObs, sizeof(double));

    dsgnmat2temptrend(tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape, tempcoeffloc,
		      tempcoeffscale, tempcoeffshape, *nSite, *nObs, usetempcov, *ntempcoeffloc,
		      *ntempcoeffscale, *ntempcoeffshape, trendlocs, trendscales, trendshapes);

    for (i=*nSite;i--;)
      for (j=*nObs;j--;)
	if (((scales[i] + trendscales[j]) <= 0) || ((shapes[i] + trendshapes[j]) <= -1)){
	  *dns = MINF;
	  return;
	}
  }

  else if (*dns != 0.0)
    return;

  //Stage 3: Transformation to unit Frechet
  if (flag)
    *dns = gev2frechTrend(data, *nObs, *nSite, locs, scales, shapes, trendlocs, trendscales,
			  trendshapes, jac, frech);

  else
    *dns = gev2frech(data, *nObs, *nSite, locs, scales, shapes, jac, frech);

  if (*dns != 0.0)
    return;
  
  if (*weighted)
    *dns = wlpliksmith(frech, mahalDist, jac, *nObs, *nSite, weights);

  else
    *dns = lpliksmith(frech, mahalDist, jac, *nObs, *nSite);

  //Stage 5: Removing the penalizing terms (if any)
  // 1- For the location parameter
  if (*locpenalty > 0)
    *dns -= penalization(locpenmat, loccoeff, *locpenalty, *nloccoeff, *npparloc);
  
  // 2- For the scale parameter
  if (*scalepenalty > 0)    
    *dns -= penalization(scalepenmat, scalecoeff, *scalepenalty, *nscalecoeff, *npparscale);
  
  // 3- For the shape parameter
  if (*shapepenalty > 0)
    *dns -= penalization(shapepenmat, shapecoeff, *shapepenalty, *nshapecoeff, *npparshape);
  
  // 4- Doing the same thing for the temporal component
  if (*temppenaltyloc > 0)
    *dns -= penalization(temppenmatloc, tempcoeffloc, *temppenaltyloc, *ntempcoeffloc,
			 *nppartempcoeffloc);

  if (*temppenaltyscale > 0)
    *dns -= penalization(temppenmatscale, tempcoeffscale, *temppenaltyscale, *ntempcoeffscale,
			 *nppartempcoeffscale);

  if (*temppenaltyshape > 0)
    *dns -= penalization(temppenmatshape, tempcoeffshape, *temppenaltyshape, *ntempcoeffshape,
			 *nppartempcoeffshape);

  // 4- Doing the same thing for the temporal component
  if (*temppenaltyloc > 0)
    *dns -= penalization(temppenmatloc, tempcoeffloc, *temppenaltyloc, *ntempcoeffloc,
			 *nppartempcoeffloc);

  if (*temppenaltyscale > 0)
    *dns -= penalization(temppenmatscale, tempcoeffscale, *temppenaltyscale, *ntempcoeffscale,
			 *nppartempcoeffscale);

  if (*temppenaltyshape > 0)
    *dns -= penalization(temppenmatshape, tempcoeffshape, *temppenaltyshape, *ntempcoeffshape,
			 *nppartempcoeffshape);

  return;
}
