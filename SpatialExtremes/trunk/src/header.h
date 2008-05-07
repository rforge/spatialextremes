#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#define RANDIN GetRNGstate()
#define RANDOUT PutRNGstate()

////////////////////////////////////
// Structure definitions
//
struct toFrech{
  int flag;
  double *frech, *jac;
};

struct covComp{
  int flag;
  double *vec;
};

struct toParam{
  int flag;
  double *locs, *scales, *shapes;
};

///////////////////////////////////
//  From schlather.c
//
void schlatherfull(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		   double *locs, double *scales, double *shapes,
		   double *scale, double *smooth, double *dns);
void schlatherlm(int *covmod, double *data, double *dist, int *nDim, int *nSite, int *nObs,
		 double *locdsgnmat, int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
		 double *shapedsgnmat, int *nshapecoeff, double *loccoeff, double *scalecoeff,
		 double *shapecoeff, double *scale, double *smooth,
		 double *dns);
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
//void smithlm(double *data, double *coord, int *nSite,
//	     int *nObs, double *locdsgnmat, int *nloccoeff,
//	     double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
//	     int *nshapecoeff, double *loccoeff, double *scalecoeff,
//	     double *shapecoeff, double *cov11, double *cov12,
//	     double *cov22, double *dns);
//void smithrb(double *data, double *distVec, int *nSite, int *degree,
//	     int *nObs, double *locdsgnmat, double *locpenmat, int *nloccoeff,
//	     double *scaledsgnmat, double *scalepenmat, int *nscalecoeff,
//	     double *shapedsgnmat, double *shapepenmat, int *nshapecoeff, 
//	     double *penalty, double *loccoeff, double *scalecoeff,
//	     double *shapecoeff, double *icov11, double *icov12,
//	     double *icov22, double *dns);
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
struct toFrech gev2frech(double *data, int nObs, int nSite, double *locs,
			 double *scales, double *shapes);
struct toParam dsgnmat2Param(double *locdsgnmat, double *scaledsgnmat,
			     double *shapedsgnmat, double *loccoeff, 
			     double *scalecoeff, double *shapecoeff,
			     int nSite, int nloccoeff, int nscalecoeff,
			     int nshapecoeff);

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
struct covComp whittleMatern(double *dist, int nPairs, double scale, double smooth);
struct covComp cauchy(double *dist, int nPairs, double scale,
		      double smooth);
struct covComp powerExp(double *dist, int nPairs, double scale,
			double smooth);
struct covComp genHyper(double *dist, int nPairs, double scale,
			double smooth1, double smooth2,
			double smooth3);
struct covComp mahalDistFct(double *distVec, int nPairs, double *icov11,
			    double *icov12, double *icov22);
struct covComp mahalDistFct2(double *distVec, int nPairs, double *cov11,
			     double *cov12, double *cov22);

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
