#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

///////////////////////////////////
//  From schlather.c
//
void schlatherfull(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		   double *locs, double *scales, double *shapes,
		   double *scale, double *smooth, double *dns);
void schlatherdsgnmat(int *covmod, double *data, double *dist, int *nDim, int *nSite, int *nObs,
		      double *locdsgnmat, double *locpenmat, int *nloccoeff, int *npparloc,
		      double *locpenalty, double *scaledsgnmat, double *scalepenmat,
		      int *nscalecoeff, int *npparscale, double *scalepenalty, double *shapedsgnmat,
		      double *shapepenmat, int *nshapecoeff, int *npparshape, double *shapepenalty,
		      double *loccoeff, double *scalecoeff, double *shapecoeff, double *scale,
		      double *smooth, double *dns);

///////////////////////////////////
//  From smith.c
//
void smithfull(double *data, double *distVec, int *nSite,
	       int *nObs, double *locs, double *scales, double *shapes,
	       double *icov11, double *icov12, double *icov22, double *dns);
void smithdsgnmat(double *data, double *distVec, int *nSite, int *nObs, 
		  double *locdsgnmat, double *locpenmat, int *nloccoeff,
		  int *npparloc, double *locpenalty, double *scaledsgnmat,
		  double *scalepenmat, int *nscalecoeff, int *npparscale,
		  double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		  int *nshapecoeff, int *npparshape, double *shapepenalty,
		  double *loccoeff, double *scalecoeff, double *shapecoeff,
		  double *icov11, double *icov12, double *icov22, double *dns);

///////////////////////////////////
//  From utils.c
//
void distance(double *coord, int *nDim, int *nSite,
	      double *dist);
void distVecFct(double *coord, int *nSite, int *nDim, double *distVec);
int gev2frech(double *data, int nObs, int nSite, double *locs,
	      double *scales, double *shapes, double *jac, double *frech);
int dsgnmat2Param(double *locdsgnmat, double *scaledsgnmat,
		  double *shapedsgnmat, double *loccoeff, 
		  double *scalecoeff, double *shapecoeff,
		  int nSite, int nloccoeff, int nscalecoeff,
		  int nshapecoeff, double *locs, double *scales,
		  double *shapes);

///////////////////////////////////
//  From univllik.c
//
void gevlik(double *data, int *n, double *loc, double *scale,
	    double *shape, double *dns);
void gpdlik(double *exceed, int *n, double *thresh, double *scale,
	    double *shape, double *dns);

///////////////////////////////////
//  From covariance.c
//
int whittleMatern(double *dist, int nPairs, double scale,
		  double smooth, double *rho);
int cauchy(double *dist, int nPairs, double scale,
	   double smooth, double *rho);
int powerExp(double *dist, int nPairs, double scale,
	     double smooth, double *rho);
int genHyper(double *dist, int nPairs, double scale,
	     double smooth1, double smooth2,
	     double smooth3, double *rho);
int mahalDistFct(double *distVec, int nPairs, double *icov11,
		 double *icov12, double *icov22, double *mahal);
int mahalDistFct2(double *distVec, int nPairs, double *cov11,
		  double *cov12, double *cov22, double *mahal);

///////////////////////////////////
//  From mcmc.c
//
SEXP gibbs(SEXP n, SEXP np, SEXP thin, SEXP init,
	   SEXP psd, SEXP f, SEXP rho);

///////////////////////////////////
//  From gradients.c
//
void smithgrad(double *data, double *distVec, int *nSite,
	       int *nObs, double *locdsgnmat, int *nloccoeff,
	       double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
	       int *nshapecoeff, double *loccoeff, double *scalecoeff,
	       double *shapecoeff, double *cov11, double *cov12,
	       double *cov22, int *fitmarge, double *grad);
void smithgrad2(double *data, double *distVec, int *nSite,
		int *nObs, double *locdsgnmat, int *nloccoeff,
		double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		int *nshapecoeff, double *loccoeff, double *scalecoeff,
		double *shapecoeff, double *icov11, double *icov12,
		double *icov22, int *fitmarge, double *grad);
void schlathergrad(int *covmod, double *data, double *dist, int *nSite,
		   int *nObs, double *locdsgnmat, int *nloccoeff,
		   double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		   int *nshapecoeff, double *loccoeff, double *scalecoeff,
		   double *shapecoeff, double *scale, double *smooth,
		   int *fitmarge, double *grad);
///////////////////////////////////
//  From pairwiselik.c
//
double lpliksmith(double *data, double *rho, double *jac,
		  int nObs, int nSite);
double lplikschlather(double *data, double *rho, double *jac,
		      int nObs, int nSite);

///////////////////////////////////
//  From penalizations.c
//
double penalization(double *penmat, double *beta, double pencoeff, int n,
		    int nppar);
double penalization2(double *penmat, double *beta, double pencoeff, int n,
		     int nppar);
