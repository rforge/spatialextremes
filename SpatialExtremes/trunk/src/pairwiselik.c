#include "header.h"

double lplikschlather(double *data, double *rho, double *jac,
		      int nObs, int nSite){
  //This function computes the log-pairwise likelihood for the
  //schalther model.

  int i, j, k, currentPair = -1;
  double c1, c2, c3, a1, a2, a3, a4, dns,
    lFvec, dvecM1, dvecM2, dvecMixed;

  for (i=0;i<(nSite-1);i++){
    for (j=i+1;j<nSite;j++){
      
      currentPair++;
      
      for (k=0;k<nObs;k++){
	c1 = 2 * data[k + i * nObs] * data[k + j * nObs] * (1 + rho[currentPair]) /
	  R_pow_di(data[k + i * nObs] + data[k + j * nObs], 2);
    
	if (c1 >= 1){
	  dns = -1.0e35;
	  return dns;
	}

	c2 = R_pow(1 - c1, .5);
    	c3 = 1/data[k + i * nObs] + 1/data[k + j * nObs];
    
	//It's the joint (log) CDF
	lFvec = exp(-c3 * (c2+1) / 2);
	
	//It's the partial derivative for marge 1
	dvecM1 = (c2+1) / 2 / R_pow_di(data[k + i * nObs], 2) -
	  c3 * (2 * c1 / (data[k + i * nObs] + data[k + j * nObs]) -
		c1 / data[k + i * nObs]) / 4 / c2;

	//if ((dvecM1 <= 0)){
	//  //printf("dvecM1 <= 0 !!!\n");
	//  dns = -1.0e35;
	//  return dns;
	//}

	//It's the partial derivative for marge 2
	dvecM2 = (c2+1) / 2 / R_pow_di(data[k + j * nObs], 2) -
	  c3 * (2 * c1 / (data[k + i * nObs] + data[k + j * nObs]) -
		c1 / data[k + j * nObs]) / 4 / c2;

	//if ((dvecM2 <= 0)){
	//  //printf("dvecM2 <= 0 !!!\n");
	//  dns = -1.0e35;
	//  return dns;
	//}
	    
	//Rmq: to have dvecM1 and dvecM2 we have to multiply
	//them by lFvec[i]. It's not done yet as dvecMixed has to be
	//computed first.

	a1 = (2 * c1 / (data[k + i * nObs] + data[k + j * nObs]) - c1 / data[k + i * nObs]) /
	  4 / R_pow_di(data[k + j * nObs], 2) / c2;
	a2 = (2 * c1 / (data[k + i * nObs] + data[k + j * nObs]) - c1 / data[k + j * nObs]) /
	  4 / R_pow_di(data[k + i * nObs], 2) / c2;
	a3 = -c3 * (-c1 / data[k + i * nObs] / data[k + j * nObs] +
		    2 * c1 / data[k + i * nObs] / (data[k + i * nObs]+data[k + j * nObs]) +
		    2 * c1 / data[k + j * nObs] / (data[k + i * nObs]+data[k + j * nObs]) -
		    6 * c1 / R_pow_di(data[k + i * nObs]+data[k + j * nObs],2)) /
	  4 / c2;
	a4 = c3 * (2 * c1 / (data[k + i * nObs]+data[k + j * nObs]) - c1 / data[k + j * nObs]) *
	  (2 * c1 / (data[k + i * nObs]+data[k + j * nObs]) - c1 / data[k + i * nObs]) / 8 /
	  R_pow(1 - c1, 1.5);
    
	//It's the mixed partial derivative
	dvecMixed = dvecM1 * dvecM2 + a1 + a2 + a3 + a4;
    
	if ((dvecMixed <= 0)){
	  //printf("dvecMixed <= 0 !!!\n");
	  dns = -1.0e35;
	  return dns;
	}

	//Now the final step, multiplying by Fvec
	//dvecM1 = log(dvecM1 * lFvec);
	//dvecM2 = log(dvecM2 * lFvec);
	dvecMixed = log(dvecMixed * lFvec) +
	  jac[k + i * nObs] + jac[k + j * nObs];

	dns = dns + dvecMixed;
      }
    }
  }

  return dns;

}

double lpliksmith(double *data, double *mahalDist, double *jac,
		  int nObs, int nSite){
  //This function computes the log-pairwise likelihood for the
  //smith model.
  
  int i, j, k, currentPair = -1;
  double c1, c2, dns, lFvec, dvecM1, dvecM2, dvecMixed;

  for (i=0;i<(nSite-1);i++){
    for (j=i+1;j<nSite;j++){
      
      currentPair++;
      
      for (k=0;k<nObs;k++){
	c1 = log(data[k + j * nObs] / data[k + i * nObs]) /
	  mahalDist[currentPair] + mahalDist[currentPair] / 2;
	c2 = log(data[k + i * nObs] / data[k + j * nObs]) /
	  mahalDist[currentPair] + mahalDist[currentPair] / 2;
	
	//It's the log of the joint CDF
	lFvec = -pnorm(c1 , 0., 1., 1, 0) / data[k + i * nObs] -
	  pnorm(c2, 0., 1., 1, 0) / data[k + j * nObs];

	//It's the partial derivative for marge 1
	dvecM1 = dnorm(c1, 0., 1., 0) / R_pow_di(data[k + i * nObs], 2) /
	  mahalDist[currentPair] - dnorm(c2, 0., 1., 0) / data[k + i * nObs] / data[k + j * nObs] /
	  mahalDist[currentPair] + pnorm(c1, 0., 1., 1, 0) / R_pow_di(data[k + i * nObs], 2);
	
	//if (dvecM1 <= 0){
	//  //printf("dvecM1 <= 0\n");
	//  dns = -1.0e35;
	//  return dns;
	//}
	
	//It's the partial derivative for marge 2
	dvecM2 = - dnorm(c1, 0., 1., 0) / data[k + i * nObs] / data[k + j * nObs] /
	  mahalDist[currentPair] + dnorm(c2, 0., 1., 0) / R_pow_di(data[k + j * nObs], 2) /
	  mahalDist[currentPair] + pnorm(c2, 0., 1., 1, 0) / R_pow_di(data[k + j * nObs], 2);
	
	//if (dvecM2 <= 0){
	//  //printf("dvecM2[%i] = %f\n", k , dvecM2[k]);
	//  dns = -1.0e35;
	//  return dns;
	//}
	
	//Rmq: to have dvecM1 and dvecM2 we have to multiply
	//them by Fvec[i]. It's not done yet as dvecMixed has to be
	//computed first.
	
	//It's the mixed partial derivative
	dvecMixed = dvecM1 * dvecM2 + c2 * dnorm(c1, 0., 1., 0) /
	  R_pow_di(data[k + i * nObs] * mahalDist[currentPair], 2) / data[k + j * nObs] +
	  c1 * dnorm(c2, 0., 1., 0) / R_pow_di(data[k + j * nObs] * mahalDist[currentPair], 2) /
	  data[k + i * nObs];
	
	if (dvecMixed <= 0){
	  //printf("dvecMixed <= 0\n");
	  dns = -1.0e35;
	  return dns;
	}
	
	//Now the final step, multiplying by Fvec
	dvecM1 = log(dvecM1) + lFvec;
	dvecM2 = log(dvecM2) + lFvec;
	dvecMixed = log(dvecMixed) + lFvec +
	  jac[k + i * nObs] + jac[k + j * nObs];

	dns = dns +  dvecMixed;
      }
    }
  }
  
  return dns;
}
