
void smithlm(double *data, double *distVec, int *nSite,
	     int *nObs, double *locdsgnmat, int *nloccoeff,
	     double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
	     int *nshapecoeff, double *loccoeff, double *scalecoeff,
	     double *shapecoeff, double *icov11, double *icov12,
	     double *icov22, double *dns){
  //This is the Smith model. It computes the pairwise log-likelihood
  
  struct toFrech toFrechObj;
  struct covComp covCompObj;
  struct toParam gevParam;
  const int nPairs = *nSite * (*nSite - 1) / 2;
  double *jac, *mahalDist, *locs, *scales, *shapes;
  
  jac = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));

  //Stage 1: Computing the Mahalanobis distance
  covCompObj = mahalDistFct(distVec, nPairs, icov11, icov12,
			    icov22);

  if (covCompObj.flag == 1){
    *dns = -1e36;
    return;
  }

  mahalDist = covCompObj.vec;

  //Stage 2: Computing the GEV parameters using the design matrix
  gevParam = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			   loccoeff, scalecoeff, shapecoeff, *nSite,
			   *nloccoeff, *nscalecoeff, *nshapecoeff);

  if (gevParam.flag == 1){
    *dns = -1e36;
    return;
  }

  locs = gevParam.locs;
  scales = gevParam.scales;
  shapes = gevParam.shapes;

  //Stage 3: Transformation to unit Frechet
  toFrechObj = gev2frech(data, *nObs, *nSite, locs, scales, shapes);
    
  if (toFrechObj.flag == 1){
    *dns = -1e36;
    return;
  }
  
  data = toFrechObj.frech;
  jac = toFrechObj.jac;
  
  //Stage 4: Bivariate density computations
  *dns = lpliksmith(data, mahalDist, jac, *nObs, *nSite);
  
}

void smithrb(double *data, double *distVec, int *nSite, int *degree,
	     int *nObs, double *locdsgnmat, double *locpenmat, int *nloccoeff,
	     double *scaledsgnmat, double *scalepenmat, int *nscalecoeff,
	     double *shapedsgnmat, double *shapepenmat, int *nshapecoeff, 
	     double *penalty, double *loccoeff, double *scalecoeff,
	     double *shapecoeff, double *icov11, double *icov12,
	     double *icov22, double *dns){
  //This is the Smith model. It computes the pairwise log-likelihood
  
  struct toFrech toFrechObj;
  struct covComp covCompObj;
  struct toParam gevParam;
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i, j;
  double *jac, *mahalDist, *locs, *scales, *shapes, pen = 0;
  
  jac = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  
  //Stage 1: Computing the Mahalanobis distance
  covCompObj = mahalDistFct(distVec, nPairs, icov11, icov12,
			    icov22);

  if (covCompObj.flag == 1){
    *dns = -1e36;
    return;
  }

  mahalDist = covCompObj.vec;

  //Stage 2: Computing the GEV parameters using the design matrix
  gevParam = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			   loccoeff, scalecoeff, shapecoeff, *nSite,
			   *nloccoeff, *nscalecoeff, *nshapecoeff);

  if (gevParam.flag == 1){
    *dns = -1e36;
    return;
  }

  locs = gevParam.locs;
  scales = gevParam.scales;
  shapes = gevParam.shapes;

  //Stage 3: Transformation to unit Frechet
  toFrechObj = gev2frech(data, *nObs, *nSite, locs, scales, shapes);
    
  if (toFrechObj.flag == 1){
    *dns = -1e36;
    return;
  }
  
  data = toFrechObj.frech;
  jac = toFrechObj.jac;
  
  //Stage 4: Bivariate density computations
  *dns = lpliksmith(data, mahalDist, jac, *nObs, *nSite);
  
  //Stage 5: Adding the penalizing terms
  // 1- For the location parameter
  for (i = 0; i < *nloccoeff; i++)
    for (j = 0; j < *nloccoeff; j++)
      pen = pen + loccoeff[j] * locpenmat[i * *nloccoeff + j] *
	loccoeff[i];

  // 2- For the scale parameter
  for (i = 0; i < *nscalecoeff; i++)
    for (j = 0; j < *nscalecoeff; j++)
      pen = pen + scalecoeff[j] * scalepenmat[i * *nscalecoeff + j] *
	scalecoeff[i];

  // 3- For the shape parameter
  for (i = 0; i < *nshapecoeff; i++)
    for (j = 0; j < *nshapecoeff; j++)
      pen = pen + shapecoeff[j] * shapepenmat[i * *nshapecoeff + j] *
	shapecoeff[i];

  *dns = *dns + R_pow_di(*penalty, *degree) * pen;

}
