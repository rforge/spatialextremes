#include "header.h"

void copula(int *copula, int *covmod, double *dist, double *data, int *nSite, int *nObs,
	    int *dim, double *locdsgnmat,  double *locpenmat, int *nloccoeff,
	    int *npparloc, double *locpenalty, double *scaledsgnmat, double *scalepenmat,
	    int *nscalecoeff, int *npparscale, double *scalepenalty,
	    double *shapedsgnmat, double *shapepenmat, int *nshapecoeff, int *npparshape,
	    double *shapepenalty, int *usetempcov, double *tempdsgnmatloc,
	    double *temppenmatloc, int *ntempcoeffloc, int *nppartempcoeffloc,
	    double *temppenaltyloc, double *tempdsgnmatscale, double *temppenmatscale,
	    int *ntempcoeffscale, int *nppartempcoeffscale, double *temppenaltyscale,
	    double *tempdsgnmatshape, double *temppenmatshape, int *ntempcoeffshape,
	    int *nppartempcoeffshape, double *temppenaltyshape, double *loccoeff,
	    double *scalecoeff, double *shapecoeff, double *tempcoeffloc,
	    double *tempcoeffscale, double *tempcoeffshape, double *DoF, double *sill, double *range,
	    double *smooth, double *smooth2, double *dns){

  int i, j, flag = usetempcov[0] + usetempcov[1] + usetempcov[2],
    nPairs = *nSite * (*nSite + 1) / 2;

  double sd = sqrt(*sill), *trendlocs, *trendscales, *trendshapes,
    *logdens = (double *) R_alloc(*nSite * *nObs, sizeof(double)),
    *covariances = (double *) R_alloc(nPairs, sizeof(double)),
    *covMat = (double *)R_alloc(*nSite * *nSite, sizeof(double)),
    *locs = (double *)R_alloc(*nSite, sizeof(double)),
    *scales = (double *)R_alloc(*nSite, sizeof(double)),
    *shapes = (double *)R_alloc(*nSite, sizeof(double)),
    *unif = (double *) R_alloc(*nSite * *nObs, sizeof(double));

  memset(covMat, 0, *nSite * *nSite * sizeof(double));

  //Stage 1: Compute the covariance at each location
  switch (*covmod){
  case 1:
    *dns = whittleMatern(dist, nPairs, *sill, *range, *smooth, covariances);
    break;
  case 2:
    *dns = cauchy(dist, nPairs, *sill, *range, *smooth, covariances);
    break;
  case 3:
    *dns = powerExp(dist, nPairs, *sill, *range, *smooth, covariances);
    break;
  case 4:
    *dns = bessel(dist, nPairs, *dim, *sill, *range, *smooth, covariances);
    break;
  case 5:
    *dns = caugen(dist, nPairs, *sill, *range, *smooth, *smooth2, covariances);
    break;
  }

  if (*dns != 0.0)
    return;

  /* We need to fill in the upper triangular part of covMatChol with
     covariances */
  {
    int current=-1;
    for (i=0;i<*nSite;i++)
      for (j=i;j<*nSite;j++){
	current++;
	covMat[i + j * *nSite] = covariances[current];
      }
  }

  //Stage 2: Compute the GEV parameters using the design matrix
  *dns = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat, loccoeff,
		       scalecoeff, shapecoeff, *nSite, *nloccoeff, *nscalecoeff,
		       *nshapecoeff, locs, scales, shapes);

  if (flag){
    trendlocs = (double *) R_alloc(*nObs, sizeof(double));
    trendscales = (double *) R_alloc(*nObs, sizeof(double));
    trendshapes = (double *) R_alloc(*nObs, sizeof(double));

    dsgnmat2temptrend(tempdsgnmatloc, tempdsgnmatscale, tempdsgnmatshape,
		      tempcoeffloc, tempcoeffscale, tempcoeffshape, *nSite,
		      *nObs, usetempcov, *ntempcoeffloc, *ntempcoeffscale,
		      *ntempcoeffshape, trendlocs, trendscales, trendshapes);

    for (i=*nSite;i--;)
      for (j=*nObs;j--;)
	if (((scales[i] + trendscales[j]) <= 0) ||
	    ((shapes[i] + trendshapes[j]) <= -1)){
	  *dns = MINF;
	  return;
	}
  }

  else if (*dns != 0.0)
    return;

  //Stage 3: Transformation to unit Frechet
  if (flag)
    *dns = gev2unifTrend(data, *nObs, *nSite, locs, scales, shapes, trendlocs,
			 trendscales, trendshapes, unif, logdens);

  else
    *dns = gev2unif(data, *nObs, *nSite, locs, scales, shapes, unif, logdens);

  if (*dns != 0.0)
    return;

  //Stage 4: Computation of the log-likelihood

  // First the "dependence" log-likelihood
  if (*copula == 1){
    for (i=(*nSite * *nObs);i--;)
      unif[i] = qnorm(unif[i], 0, sd, 1, 0);

    *dns = gaussianCopula(unif, sd, covMat, *nObs, *nSite);
  }

  if (*copula == 2){

    if (*DoF <= 0){
      *dns = (1 - *DoF) * (1 - *DoF) * MINF;
    }
      
    for (i=(*nSite * *nObs);i--;)
      unif[i] = qt(unif[i], *DoF, 1, 0);

    *dns = studentCopula(unif, *DoF, covMat, *nObs, *nSite);
  }

  if (*dns == MINF)
    return;

  /* Second the "marginal" log-likelihood i.e., log of the GEV
  densities */
  for (i=(*nSite * *nObs);i--;)
    *dns += logdens[i];

  //Stage 5: Removing the penalizing terms (if any)
  // 1- For the location parameter
  if (*locpenalty > 0)
    *dns -= penalization(locpenmat, loccoeff, *locpenalty, *nloccoeff,
			 *npparloc);
  
  // 2- For the scale parameter
  if (*scalepenalty > 0)    
    *dns -= penalization(scalepenmat, scalecoeff, *scalepenalty, *nscalecoeff,
			 *npparscale);
  
  // 3- For the shape parameter
  if (*shapepenalty > 0)
    *dns -= penalization(shapepenmat, shapecoeff, *shapepenalty, *nshapecoeff,
			 *npparshape);

  // 4- Doing the same thing for the temporal component
  if (*temppenaltyloc > 0)
    *dns -= penalization(temppenmatloc, tempcoeffloc, *temppenaltyloc,
			 *ntempcoeffloc, *nppartempcoeffloc);

  if (*temppenaltyscale > 0)
    *dns -= penalization(temppenmatscale, tempcoeffscale, *temppenaltyscale,
			 *ntempcoeffscale, *nppartempcoeffscale);

  if (*temppenaltyshape > 0)
    *dns -= penalization(temppenmatshape, tempcoeffshape, *temppenaltyshape,
			 *ntempcoeffshape, *nppartempcoeffshape);

  if (!R_FINITE(*dns))
    *dns = MINF;

  return;
  
}

double gaussianCopula(double *data, double sd, double *covMat, int nObs,
		      int nSite){
  // This function computes the log-likelihood for the Gaussian copula 

  int i, j, info = 0, oneInt = 1;
  double ans = 0, logDet = 0, one = 1;

  // Cholesky decomposition
  F77_CALL(dpotrf)("U", &nSite, covMat, &nSite, &info);

  if (info != 0)
    return MINF;
  
  // Compute the log of the determinant from this Cholesky decomp
  for (i=nSite;i--;)
    logDet += log(covMat[i * (nSite + 1)]);
  
  logDet *= 2;

  ans = - 0.5 * nObs * (logDet + nSite * log(M_2PI));

  // Compupte the log of the multivariate Gaussian density
  double *dummy = (double *) R_alloc(nSite, sizeof(double));
  for (i=nObs;i--;){
    for (j=nSite;j--;)
      dummy[j] = data[i + j * nObs];

    F77_CALL(dtrsm)("L", "U", "T", "N", &nSite, &oneInt, &one, covMat,
		    &nSite, dummy, &nSite);

    for (j=nSite;j--;)
      ans -= 0.5 * dummy[j] * dummy[j];
  }

  //Jacobian part of the gaussian copula
  for (i=(nSite * nObs);i--;)
    ans -= dnorm(data[i], 0, sd, 1);

  return ans;
}


double studentCopula(double *data, double DoF, double *covMat, int nObs,
		     int nSite){
  // This function computes the log-likelihood for the Student copula 

  int i, j, info = 0, oneInt = 1;
  double ans = 0, logDet = 0, one = 1;

  // Cholesky decomposition
  F77_CALL(dpotrf)("U", &nSite, covMat, &nSite, &info);

  if (info != 0)
    return MINF;
  
  // Compute the log of the determinant from this Cholesky decomp
  for (i=nSite;i--;)
    logDet += log(covMat[i * (nSite + 1)]);
  
  logDet *= 2;

  // Compupte the log of the multivariate Student density
  double *dummy = (double *) R_alloc(nSite, sizeof(double));
  for (i=nObs;i--;){
    for (j=nSite;j--;)
      dummy[j] = data[i + j * nObs];

    F77_CALL(dtrsm)("L", "U", "T", "N", &nSite, &oneInt, &one, covMat,
		    &nSite, dummy, &nSite);

    for (j=nSite;j--;)
      ans += dummy[j] * dummy[j];
  }

  ans = nObs * (lgammafn(0.5 * (DoF + nSite)) - lgammafn(0.5 * DoF) - 0.5 *
		nSite * log(M_PI * DoF) - 0.5 * logDet) -
    0.5 * (DoF + nSite) * log1p(ans / DoF);

  //Jacobian part of the Student copula
  for (i=(nSite * nObs);i--;)
    ans -= dt(data[i], DoF, 1);

  return ans;
}

