#include "header.h"

void rbrowndirect(double *coord, double *bounds, int *nObs, int *nSite,
		  int *dim, int *grid, double *range, double *smooth,
		  double *ans){
  /* This function generates random fields for the geometric model

     coord: the coordinates of the locations
      nObs: the number of observations to be generated
     nSite: the number of locations
       dim: the random field is generated in R^dim
      grid: Does coord specifies a grid?
     range: the range parameter
    smooth: the smooth parameter
       ans: the generated random field */

  int neffSite, covmod = 6, lagi = 1, lagj = 1, nSimTot = 250, oneInt = 1;
  const double irange = 1 / *range, inflate = R_pow(4, 1 / *smooth) * *range;
  double one = 1, zero = 0;

  //Inflate the region
  for (int i=0;i<*dim;i++){
    bounds[2 * i] -= inflate;
    bounds[2 * i + 1] += inflate;
  }

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  GetRNGstate();
  for (int i=0; i<*nObs; i++){
    double poisson = 0;
    int nKO = neffSite;
    
    double *gp = (double *)R_alloc(neffSite, sizeof(double)),
      *covmat = (double *)R_alloc(neffSite * neffSite, sizeof(double));

    while (nKO){
    //for (int nSim=nSimTot; nSim--;){
      double *shift = (double *) R_alloc(*dim, sizeof(double));
      double *shiftedCoord = (double *) R_alloc(*dim * *nSite, sizeof(double));
      double *vario = (double *) R_alloc(neffSite, sizeof(double));

      // Shift the locations
      for (int j=0; j<*dim;j++)
	shift[j] = runif(bounds[2 * j], bounds[2 * j + 1]);
      
      for (int j=0; j<*nSite;j++)
	for (int k=0; k<*dim; k++)
	  shiftedCoord[k * *nSite + j] = coord[k * *nSite + j] - shift[k];
      
      // Compute the variogram
      distance2orig(shiftedCoord, *nSite, *dim, vario, *grid);
      
      for (int j=0; j<neffSite; j++)
	vario[j] = R_pow(vario[j] * irange, *smooth);
      
      // Compute the covariance matrix for the shifted locations
      buildcovmat(nSite, grid, &covmod, shiftedCoord, dim, &zero, &one, range,
		  smooth, covmat);
      
      /* Compute the Cholesky decomposition of the covariance matrix. */
      int info = 0;
      F77_CALL(dpotrf)("U", &neffSite, covmat, &neffSite, &info);

      if (info != 0)
	error("error code %d from Lapack routine '%s'", info, "dpotrf");

      /* Simulate a std normal random vector from this decomposition */
      for (int j=0; j<neffSite;j++)
	gp[j] = norm_rand();

      F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite, gp, &oneInt);
      
      poisson += exp_rand();
      double ipoisson = -log(poisson), thresh = 5 + ipoisson;
	
      nKO = neffSite;
      for (int j=0; j<neffSite; j++){
	ans[j * lagj + i * lagi] = fmax2(gp[j] - vario[j] + ipoisson,
					 ans[j * lagj + i * lagi]);
	nKO -= (thresh <= ans[j * lagj +  i * lagi]);
      }
    }
  }

  PutRNGstate();

  for (int i=0; i<(*nObs * neffSite);i++)
    ans[i] = exp(ans[i]);

  return;
}
