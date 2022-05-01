#ifndef TYPES_AND_CONSTANTS_H
#define TYPES_AND_CONSTANTS_H

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multiroots.h>

#define N_POINTS 100 //number of points

typedef struct potentialParams{
	int type;
	double v0;
	double a; // first point after 0, it's the upper bound of the 1st domain (where V=0);
	double b; // second point, it's the interface between the 2nd and 3rd domains (respectively where V=V0 and V=0)
}potentialParams;

typedef struct schrodingerParameters{
	potentialParams potential;
	double energy;
	double bound; // third point, it's the upper bound of the 3d domain (after that point V=infinity)
	double alpha;
	int doDraw; // is equal to 0 if we don't want to draw phi(x) on the graph
	int currentDomain; // used to know in which domain we are to solve the ODE
	double prevDomainY0; // phi(x) at the right bound of the previous domain
	double prevDomainY1; // dphi(x)/dx at the right bound of the previous domain
	double prevDomainY2; // N(x) at the right bound of the previous domain
}schrodingerParameters;

#endif