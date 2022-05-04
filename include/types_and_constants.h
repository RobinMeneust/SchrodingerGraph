/**
 * \file types_and_constants.h
 * \brief Defines the custom types and constants used
 * \date 2022
 */

#ifndef TYPES_AND_CONSTANTS_H
#define TYPES_AND_CONSTANTS_H

/**
 * \def N_POINTS 
 * \brief Constant corresponding to the number of number of points drawn in the graph for phi(x)
 */
#define N_POINTS 500

/**
 * \struct potentialParams
 * \brief Contains data about the potential. Those parameters are sent to the solver
 */

typedef struct potentialParams{
	int type; /*!< type of the potential (e.g 0 if it's null) */
	double v0; /*!< value of the potential in the well if it isn't null */
	double a; /*!< first point after 0 in the graph, it's the upper bound of the 1st domain (where V=0); */
	double b; /*!< second point in the graph, it's the interface between the 2nd and 3rd domains (respectively where V=V0 and V=0) */
}potentialParams;

/**
 * \struct schrodingerParameters
 * \brief Contains all the parameters necessary for the solver to obtain a solution to the equation
 */

typedef struct schrodingerParameters{
	potentialParams potential; /*!< contains data about the potential */
	double energy; /*!< energy of the particule in the well */
	double bound; /*!< third point in the graph, it's the upper bound of the 3d domain (after that point V=infinity) */
	double alpha; /*!< the coefficient for phi(x) in the equation: phi''(x) + alpha * phi(x) = 0 (i.e alpha = 2*m / (hbar*hbar)) */
	int doDraw; /*!< is equal to 0 if we don't want to draw phi(x) on the graph */
	int currentDomain; /*!< used to know in which domain we are going to solve the ODE */
	double prevDomainY0; /*!< phi(x) at the right bound of the previous domain */
	double prevDomainY1; /*!< dphi(x)/dx at the right bound of the previous domain */
	double prevDomainY2; /*!< N(x) at the right bound of the previous domain */
}schrodingerParameters;

#endif