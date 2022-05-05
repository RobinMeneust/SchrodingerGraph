/**
 * \file schrodinger_functions.h
 * \brief Contains the functions prototypes and the required includes of schrodinger_functions.c
 * \date 2022
 */

#ifndef SCHRODINGER_FUNCTIONS
#define SCHRODINGER_FUNCTIONS

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multiroots.h>
#include "../include/types_and_constants.h"

double getPotential(double x, potentialParams potential);
double getCoefficientFromParams(double x, schrodingerParameters params);
int y_derivative (double x, const double y[], double f[], void *params);
int solveODE(double z, schrodingerParameters params, double f[3]);
double findRoot(schrodingerParameters params);
int solveODEMultipleDomains(const gsl_vector* input, void* params, gsl_vector* f);
int findMultipleRoots(double x_init, schrodingerParameters params, double roots[2]);
void savePotential(schrodingerParameters params);
void solveSchrodinger(schrodingerParameters* params);

#endif