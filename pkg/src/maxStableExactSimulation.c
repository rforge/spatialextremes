#include "header.h"

void rschlatherexact(double *coord, int *nObs, int *nSite, int *dim,
		     int *covmod, int *grid, double *nugget, double *range,
		     double *smooth, double *ans){

  /*
     This function generates random fields for the Schlather model
     using the exact procedure of Dombry et al. (2106) "Exact
     simulation of max-stable processes" Biometrika

     coord: the coordinates of the locations
      nObs: the number of observations to be generated
    nSite: the number of locations
       dim: the random field is generated in R^dim
    covmod: the covariance model
      grid: Does coord specifies a grid?
      nugget: the nugget parameter
     range: the range parameter
    smooth: the smooth parameter
       ans: the generated random field
  */

  int neffSite, lagi = 1, lagj = 1, oneInt = 1;
  double sill = 1 - *nugget;
  const double normCst = M_SQRT2 * M_SQRT_PI;


  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  double *covmat = malloc(neffSite * neffSite * sizeof(double)),
    *gp = malloc(neffSite * sizeof(double));

  buildcovmat(nSite, grid, covmod, coord, dim, nugget, &sill, range,
	      smooth, covmat);

  /* Compute the Cholesky decomposition of the covariance matrix once for all */
  int info = 0;
  F77_CALL(dpotrf)("U", &neffSite, covmat, &neffSite, &info);

  if (info != 0)
    error("error code %d from Lapack routine '%s'", info, "dpotrf");

  GetRNGstate();

  for (int i=0; i<*nObs; i++){

    // 1. Generate the first extremal function
    for (int j=0;j<neffSite;j++)
      gp[j] = norm_rand();

    F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite, gp, &oneInt);
    for (int j=0;j<neffSite;j++)
      gp[j] *= sign(gp[0]) * normCst;//sign to ensure that gp[0] > 0

    for (int j=0;j<neffSite; j++)
      ans[j * lagj + i * lagi] = fmax2(gp[j] / exp_rand(), 0.0);

    // Rmk: 0.005 is a numerical trick to avoid infinite looping

    // 2. Generate the remaining extremal functions (in a iterative way)
    for (int j=1;j<neffSite;j++){
      double poisson = 0, ipoisson = R_PosInf;

      do {
	// Simulate the Poisson process on (0, \infty)
	poisson += exp_rand();
	ipoisson = 1 / poisson;

	// Simulate the Y process (normalized such that Y(s_j) = 1)
	for (int k=0;k<neffSite;k++)
	  gp[k] = norm_rand();

	F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite, gp, &oneInt);

	for (int k=0;k<neffSite; k++)
	  gp[k] *= sign(gp[j]) * normCst;//sign to ensure that gp[j] > 0


	for (int k=0;k<j;k++){
	  if ((ipoisson * gp[k]) > ans[k * lagj + i * lagi])
	    break;

	  if (k == (j-1)){
	    /*
	      If true then it is a valid extremal function (since it
	      does not exceed the previous maxima at location coord_1,
	      ..., coord_{j-1}
	    */

	    for (int l=j;l<neffSite;l++)
	      ans[l * lagj + i * lagi] = fmax2(ans[l * lagj + i * lagi], ipoisson * gp[l]);
	  }
	}
      } while ((ipoisson * normCst * 3.5) > ans[j * lagj + i * lagi]);
    }
  }

  PutRNGstate();
  free(covmat); free(gp);

  return;
}

void rbrownrexact(double *coord, int *nObs, int *nSite, int *dim,
		  int *covmod, int *grid, double *range, double *smooth,
		  double *ans){

  /*
     This function generates random fields for the Brown-Resnick model
     using the exact procedure of Dombry et al. (2106) "Exact
     simulation of max-stable processes" Biometrika

     coord: the coordinates of the locations
      nObs: the number of observations to be generated
    nSite: the number of locations
       dim: the random field is generated in R^dim
    covmod: the covariance model
      grid: Does coord specifies a grid?
     range: the range parameter
    smooth: the smooth parameter
       ans: the generated random field
  */

  int neffSite, lagi = 1, lagj = 1, oneInt = 1, zeroInt = 0;
  double zero = 0, one = 1, irange = 1 / *range;
  *covmod = 6;

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  double *covmat = malloc(neffSite * neffSite * sizeof(double)),
    *gp = malloc(neffSite * sizeof(double)),
    *vario = malloc(neffSite * sizeof(double)),
    *shiftedCoord = malloc(*nSite * *dim * sizeof(double)),
    *orig = malloc(*dim * sizeof(double));

  buildcovmat(nSite, grid, covmod, coord, dim, &zero, &one, range,
	      smooth, covmat);

  /* Compute the Cholesky decomposition of the covariance matrix once for all */
  int info = 0;
  F77_CALL(dpotrf)("U", &neffSite, covmat, &neffSite, &info);

  if (info != 0)
    error("error code %d from Lapack routine '%s'", info, "dpotrf");

  GetRNGstate();

  for (int i=0; i<*nObs; i++){

    for (int d=0;d<*dim;d++)
      orig[d] = coord[d * *nSite];

    // 1. Generate the first extremal function
    for (int j=0;j<neffSite;j++)
      gp[j] = norm_rand();

    F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite, gp, &oneInt);

    // Compute the variogram where the origin is fixed at coord[0]
    for (int j=0; j<*nSite;j++)
      for (int k=0; k<*dim; k++)
	shiftedCoord[k * *nSite + j] = coord[k * *nSite + j] - orig[k];

    distance2orig(shiftedCoord, *nSite, *dim, vario, *grid);

    for (int j=0; j<neffSite; j++)
      vario[j] = R_pow(vario[j] * irange, *smooth);

    for (int j=0;j<neffSite; j++)
      ans[j * lagj + i * lagi] = exp(gp[j] - gp[0] - vario[j]) / exp_rand();

    // 2. Generate the remaining extremal functions (in a iterative way)
    for (int j=1;j<neffSite;j++){
      double poisson = 0, ipoisson = R_PosInf;

      // Set the origin to the j-th location
      for (int d=0;d<*dim;d++)
	orig[0] = coord[j / *nSite];
      orig[1] = coord[j - (j / *nSite) * *nSite];

      do {
	// Simulate the Poisson process on (0, \infty)
	poisson += exp_rand();
	ipoisson = 1 / poisson;

	// Simulate the Y process (normalized such that Y(s_j) = 1)
	for (int k=0;k<neffSite;k++)
	  gp[k] = norm_rand();

	F77_CALL(dtrmv)("U", "T", "N", &neffSite, covmat, &neffSite, gp, &oneInt);


	// Compute the variogram where the origin is fixed at coord[0]
	for (int k=0; k<*nSite;k++)
	  for (int l=0; l<*dim; l++)
	    shiftedCoord[l * *nSite + k] = coord[l * *nSite + k] - orig[l];

	distance2orig(shiftedCoord, *nSite, *dim, vario, *grid);

	for (int k=0; k<neffSite; k++)
	  vario[k] = R_pow(vario[k] * irange, *smooth);

	for (int k=0;k<neffSite; k++)
	  gp[k] = exp(gp[k] - gp[j] - vario[k]);


	for (int k=0;k<j;k++){
	  if ((ipoisson * gp[k]) > ans[k * lagj + i * lagi])
	    break;

	  if (k == (j-1)){
	    /*
	      If true then it is a valid extremal function (since it
	      does not exceed the previous maxima at location coord_1,
	      ..., coord_{j-1}
	    */

	    for (int l=j;l<neffSite;l++)
	      ans[l * lagj + i * lagi] = fmax2(ans[l * lagj + i * lagi], ipoisson * gp[l]);
	  }
	}
      } while (ipoisson > ans[j * lagj + i * lagi]);
    }
  }

  PutRNGstate();
  free(covmat); free(gp); free(vario); free(shiftedCoord); free(orig);

  return;
}
