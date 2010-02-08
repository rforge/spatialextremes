#include "header.h"

void circemb(int *nsim, int *ngrid, double *steps, int *dim, int *covmod,
	     double *nugget, double *sill, double *range, double *smooth,
	     double *ans){

  int m = 2, i, j, k, r, nbar = *ngrid * *ngrid;
  //irho is the imaginary part of the covariance -> 0
  double *rho, *irho; 
    
  /* Find the smallest size m for the circulant embedding matrix */
  {
    int dummy = 2 * (*ngrid - 1);
    do {
      m *= 2;
    } while (m < dummy);
  }
  
  /* ---------- beginning of the embedding stage ---------- */
  int mbar = m * m, halfM = m / 2, notPosDef = 0;
  do {
    double *dist;
    dist = (double *)R_alloc(mbar, sizeof(double));

    notPosDef = 0;
    //Computation of the distance
    for (r=mbar;r--;){
      i = r % m;
      j = r / m;
      
      if (i > halfM)
	i -= m;
      
      if (j > halfM)
	j -= m;
      
      dist[r] = pythag(steps[0] * i, steps[1] * j);
    }

    //Computations of the covariances
    rho = (double *)R_alloc(mbar, sizeof(double));
    irho = (double *)R_alloc(mbar, sizeof(double));
    memset(irho, 0, mbar * sizeof(double));

    switch (*covmod){
    case 1:
      whittleMatern(dist, mbar, *sill, *range, *smooth, rho);
      break;
    case 2:
      cauchy(dist, mbar, *sill, *range, *smooth, rho);
      break;
    case 3:
      powerExp(dist, mbar, *sill, *range, *smooth, rho);
      break;
    case 4:
      bessel(dist, mbar, *dim, *sill, *range, *smooth, rho);
      break;
    }

    /* Compute the eigen values to check if the circulant embbeding
       matrix is positive definite */

    /* Note : The next lines is only valid for 2d random fields. I
       need to change if there are m_1 \neq m_2 as I suppose that m_1
       = m_2 = m */
    int maxf, maxp, *iwork;
    double *work;

    fft_factor(m, &maxf, &maxp);
    work = (double *)R_alloc(4 * maxf, sizeof(double));
    iwork = (int *)R_alloc(maxp, sizeof(int));
    fft_work(rho, irho, m, m, 1, -1, work, iwork);

    fft_factor(m, &maxf, &maxp);
    work = (double *)R_alloc(4 * maxf, sizeof(double));
    iwork = (int *)R_alloc(maxp, sizeof(int));
    fft_work(rho, irho, 1, m, m, -1, work, iwork);

    //Check if the eigenvalues are all positive
    for (i=mbar;i--;){
      notPosDef |= (rho[i] <= 0) || (fabs(irho[i]) > 0.001);
    }

    if (notPosDef){
      halfM = m;
      m *= 2;
      mbar = m * m;
    }

    if (mbar > 1048576)
      error("Impossible to embbed the covariance matrix");
    
  } while (notPosDef);
  /* --------- end of the embedding stage --------- */

  /* Computation of the square root of the eigenvalues */
  for (i=mbar;i--;){
    rho[i] = sqrt(rho[i]);
    irho[i] = 0;//No imaginary part
  }

  /* The method generates two independent random fields per iteration
     so there's no need to perform nsim iterations but only
     neffSim. */
  int neffSim = *nsim / 2 + *nsim % 2, mdag = m / 2 + 1,
    mdagbar = mdag * mdag;
  double *a, *ia, factor = 1 / sqrt(mbar);

  a = (double *)R_alloc(mbar, sizeof(double));
  ia = (double *)R_alloc(mbar, sizeof(double));
    
  GetRNGstate();
  for (k=*nsim;k--;){
    
    /* ---------- Simulation from \Lambda^1/2 Q* Z ------------ */
    for (r=mdagbar;r--;){
      /* Below is the procedure 5.2.4 in Wood and Chan */

      //Computation of the cardinality of A(j)
      int i = r % mdag, j = r / mdag, j1, j2;
      double u, v;

      int card = (i != 0) * (i != halfM) + 2 * (j != 0) * (j != halfM);
      
      switch (card){
      case 0:
	j1 = i + m * j;
	a[j1] = factor * rho[j1] * norm_rand();
	ia[j1] = 0;
	break;
      case 1:
	//B(1) = 0, B^c(1) = {1}
	factor *= M_SQRT1_2;
	j1 = i + m * j;
	j2 = m - i + m * j;
	u = norm_rand();
	v = norm_rand();
	a[j1] = ia[j1] = factor * rho[j1];
	a[j1] *= u;
	ia[j1] *= v;
	a[j2] = ia[j2] = factor * rho[j2];
	a[j2] *= u;
	ia[j2] *= -v;
	break;
      case 2:
	//B(1) = 0, B^c(1) = {2}
	factor *= M_SQRT1_2;
	j1 = i + m * j;
	j2 = i + m * (m - j);
	u = norm_rand();
	v = norm_rand();
	a[j1] = ia[j1] = factor * rho[j1];
	a[j1] *= u;
	ia[j1] *= v;
	a[j2] = ia[j2] = factor * rho[j2];
	a[j2] *= u;
	ia[j2] *= -v;
	break;
      case 3:
	//B(1) = {1}, B^c(1) = {2}
	factor *= M_SQRT1_2;
	j1 = (m - i) + m * j;
	j2 = i + m * (m - j);
	u = norm_rand();
	v = norm_rand();
	a[j1] = ia[j1] = factor * rho[j1];
	a[j1] *= u;
	ia[j1] *= v;
	a[j2] = ia[j2] = factor * rho[j2];
	a[j2] *= u;
	ia[j2] *= -v;
	
	//B(2) = {1,2}, B^c(2) = {0}
	j1 = (m - i) + m * (m - j);
	j2 = i + m * j;
	u = norm_rand();
	v = norm_rand();
	a[j1] = ia[j1] = factor * rho[j1];
	a[j1]*= u;
	ia[j1] *= v;
	a[j2] = ia[j2] = factor * rho[j2];
	a[j2]*= u;
	ia[j2] *= -v;      
	break;
      }
    }

    /* ---------- Computation of Q \Lambda^1/2 Q* Z ------------ */
    int maxf, maxp, *iwork;
    double *work;
    
    /* The next lines is only valid for 2d random fields. I need to
       change if m_1 \neq m_2 as here I suppose that m_1 = m_2 = m */
    fft_factor(m, &maxf, &maxp);
    work = (double *)R_alloc(4 * maxf, sizeof(double));
    iwork = (int *)R_alloc(maxp, sizeof(int));
    fft_work(a, ia, m, m, 1, -1, work, iwork);
    
    fft_factor(m, &maxf, &maxp);
    work = (double *)R_alloc(4 * maxf, sizeof(double));
    iwork = (int *)R_alloc(maxp, sizeof(int));
    fft_work(a, ia, 1, m, m, -1, work, iwork);
    
    /* ------------ Get the random field ------------ */
    /*if ((2 * (k + 1)) < *nsim){ 
      int idx1, idx2;
      for (i=nbar;i--;){
	idx1 = i % m;
	idx2 = i / m;
	ans[i + 2 * k * nbar] = a[idx1 + m * idx2];
	ans[i + (2 * k + 1) * nbar] = ia[idx1 + m * idx2];
      }
    }

    else 
      for (i=nbar;i--;)
	ans[i + 2 * k * nbar] = a[(i % m) + m * (i / m)];
    */

    for (i=nbar;i--;)
      ans[i + k * nbar] = a[i % *ngrid + m * (i / *ngrid)];
  }
  PutRNGstate();  
    
  if (*nugget > 0){
    int dummy = *nsim * nbar;
    double sqrtNugget = sqrt(*nugget);
    
    GetRNGstate();
    for (i=dummy;i--;)
      ans[i] += sqrtNugget * norm_rand();

    PutRNGstate();
  }

  return;
}
