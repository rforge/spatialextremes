#include "header.h"

void rsmith(double *coord, double *center, double *edge, int *nObs,
	    int *nSites, int *grid, double *cov11, double *cov12,
	    double *cov22, double *ans){
  /* coord: the coordinates of the locations
    center: the center of the compact set - here I use a square
      edge: the length of the edge of the square
      nObs: the number of observations to be generated
      grid: Does coord specifies a grid?
    nSites: the number of locations
     covXX: the parameters of the bivariate normal density
       ans: the generated random field */

  const double det = *cov11 * *cov22 - *cov12 * *cov12,
    uBound = 1 / (M_2PI * sqrt(det)), itwiceDet = 1 / (2 * det);
  int i, j, k, nKO, inSites;
  double poisson, ipoisson, u1, u2, y, theta, r, thresh, lebesgue;

  if ((det <= 0) || (*cov11 <= 0))
    error("The covariance matrix isn't semi-definite positive!\n");

  /* We first center the coordinates to avoid repeation of unnecessary
    operations in the while loop */
  for (i=0;i<*nSites;i++){
    coord[i] -= center[0];
    coord[*nSites + i] -= center[1];
  }

  /* Simulation according to the Schlather methodology.  The compact
     set need to be inflated first */
  *edge += 3.46 * sqrt(fmax2(*cov11, *cov22));
  lebesgue = *edge * *edge;

  GetRNGstate();
  
  if (*grid){
    //Simulation part if a grid is used
    for (i=*nObs;i--;){
      inSites = i * *nSites * *nSites;
      poisson = 0;
      nKO = *nSites * *nSites;
    
      while (nKO) {
	/* The stopping rule is reached when nKO = 0 i.e. when each site
	   satisfies the condition in Eq. (8) of Schlather (2002) */

	poisson += exp_rand();
	ipoisson = 1 / poisson;
	thresh = uBound * ipoisson;

	//We simulate points uniformly in [-r/2, r/2]^2
	u1 = *edge * runif(-0.5, 0.5);
	u2 = *edge * runif(-0.5, 0.5);
      
	nKO = *nSites * *nSites;
	for (j=*nSites;j--;){
	  for (k=*nSites;k--;){
	    if (thresh > ans[j + k * *nSites + inSites]){
	      /* This is the bivariate normal density with 0 mean and
		 cov. matrix [cov11, cov12; cov12, cov22] */
	      y = exp((-*cov22 * (coord[j] - u1) * (coord[j] - u1) + 2 * *cov12 *
		       (coord[j] - u1) * (coord[*nSites + k] - u2) - *cov11 *
		       (coord[*nSites + k] - u2) * (coord[*nSites + k] - u2)) *
		      itwiceDet) * thresh;
	      
	      ans[j + k * * nSites + inSites] = fmax2(y, ans[j + k * *nSites + inSites]);
	    }
	
	    else
	      nKO--;
	  }
	}
      }
    }
  }

  else{
    //Simulation part if a grid isn't used
    for (i=*nObs;i--;){
      inSites = i * *nSites;
      poisson = 0;
      nKO = *nSites;
    
      while (nKO) {
	/* The stopping rule is reached when nKO = 0 i.e. when each site
	   satisfies the condition in Eq. (8) of Schlather (2002) */

	poisson += exp_rand();
	ipoisson = 1 / poisson;
	thresh = uBound * ipoisson;

	//We simulate points uniformly in [-r/2, r/2]^2
	u1 = *edge * runif(-0.5, 0.5);
	u2 = *edge * runif(-0.5, 0.5);
      
	nKO = *nSites;
	for (j=*nSites;j--;){
	  if (thresh > ans[j + inSites]){
	    /* This is the bivariate normal density with 0 mean and
	       cov. matrix [cov11, cov12; cov12, cov22] */
	    y = exp((-*cov22 * (coord[j] - u1) * (coord[j] - u1) + 2 * *cov12 *
		     (coord[j] - u1) * (coord[*nSites + j] - u2) - *cov11 *
		     (coord[*nSites + j] - u2) * (coord[*nSites + j] - u2)) *
		    itwiceDet) * thresh;
	  
	    ans[j + inSites] = fmax2(y, ans[j + inSites]);
	  }
	
	  else
	    nKO--;
	}
      }
    }
  }
 
  PutRNGstate();

  /* Lastly, we multiply by the Lebesgue measure of the dilated
    compact set */
  if (*grid){
    for (i=(*nSites * *nSites * *nObs);i--;)
      ans[i] *= lebesgue;
  }
  
  else{
    for (i=(*nSites * *nObs);i--;)
      ans[i] *= lebesgue;
  }

  return;
}
