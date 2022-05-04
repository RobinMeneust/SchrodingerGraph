#include "../include/types_and_constants.h"
#include "../include/schrodinger_functions.h"

double getPotential(double x, potentialParams potential){
	switch(potential.type){
		case 0: 	// V=0
			return 0;
		case 1:		// step potential
			if(x<potential.a){
				return 0;
			}
			else{
				return potential.v0;
			}
		case 2:		// rectangular potential
			if(x<potential.a || x>potential.b){
				return 0;
			}
			else{
				return potential.v0;
			}
		default:
			fprintf(stderr, "ERROR : in getPotential(), parameters are incorrect\n");
			exit(EXIT_FAILURE);
	}
}


// Calculate the coefficient of phi (C) which is in the equation of Schrodinger : phi'' + C*phi = 0
double getCoefficientFromParams(double x, schrodingerParameters params)
{
	return params.alpha*(params.energy-getPotential(x, params.potential));
}

// calculate f(y0,y1,y2)
int y_derivative (double x, const double y[], double f[] /* = dydt*/, void *params){
	double alpha = getCoefficientFromParams(x, *(schrodingerParameters*) params);
	f[0]=y[1];
	f[1]=(-1.0)*alpha*y[0];
	f[2]=y[0]*y[0];

	return GSL_SUCCESS;
}

// solve the ODE with the given parameters to get the function phi(x), and draw it if params.doDraw!=0

int solveODE(double z, schrodingerParameters params, double f[3]){
	FILE* dataFile=NULL;

	double x=0.0;
	double length=params.bound;
	double y1_ini=0.0, y2_ini=z, y3_ini=0.0; // phi(0)=0 and N(0)=0

	if(params.doDraw)
		dataFile = fopen("data/phi.dat", "a");
	
	gsl_odeiv2_system sys = {y_derivative, NULL, 3, &params}; // we initialize the ODE system
	// we set the bounds
	if(params.potential.type!=0 && params.currentDomain==0){
		x = 0.0;
		length=params.potential.a-x;
	}
	else if(params.currentDomain==1){
		x = params.potential.a;
		length=params.potential.b-x;
		y1_ini=params.prevDomainY0;
		y2_ini=params.prevDomainY1;
		y3_ini=params.prevDomainY2;
	}
	else if(params.currentDomain==2){
		x = params.potential.b;
		length=params.bound-x;
		y1_ini=params.prevDomainY0;
		y2_ini=params.prevDomainY1;
		y3_ini=params.prevDomainY2;
	}


	gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

	double y[3] = {y1_ini, y2_ini, y3_ini}; //(y0,y1,y2)(x=0) = (0,z,0)
	if(params.doDraw){
		fprintf(dataFile, "%lf %lf\n", x, y[0]);
	}
	//We get the values of (y0,y1,y2)
	double xi=x;

	for (int i = 1; i <= N_POINTS; i++)
	{
		xi += length / (double)N_POINTS;
		//printf("xi : %lf\n", xi);
		//if(xi>0.35)
		//	sleep(1);
		int status = gsl_odeiv2_driver_apply(driver, &x, xi, y);
		if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}
		//printf("f(%.2lf)=%lf\n", x, y[0]);
		//printf("Z=%lf | Y(%lf) = %lf %lf %lf\n", z, x, y[0], y[1], y[2]);
		if(params.doDraw)
			fprintf(dataFile, "%lf %lf\n", x, y[0]);
	}

	if(f!=NULL){
		f[0]=y[0];
		f[1]=y[1];
		f[2]=y[2];
	}

	gsl_odeiv2_driver_free (driver);
	if(params.doDraw)
		fclose(dataFile);
	
	return GSL_SUCCESS;
}

// Find the root with the shooting method
double findRoot(schrodingerParameters params){
	double y[3]={0.0, 0.0, 0.0};
	double step=0.01;
	double z=0.0-step;
	double eps=1.0;
	int iter=0;
	double best_z_approx=z;
	double eps_min=eps;

	do{
		z+=step;
		solveODE(z, params, y);
		eps=fabs(y[2]-1); //N(L)-1
		if(eps<eps_min){
			eps_min=eps;
			best_z_approx=z;
		}
		iter++;
	}while(eps>1e-8 && iter<5000);
	
	return best_z_approx;
}

// call solveODE() to solve the equation on different domains with different potentials
int solveODEMultipleDomains(const gsl_vector* input, void* params, gsl_vector* f){
	schrodingerParameters parameters = *(schrodingerParameters*) params;
	double y1[3]={0.0};
	double y2[3]={0.0};
	double y3[3]={0.0};
	double z=gsl_vector_get(input, 1);
	parameters.energy=gsl_vector_get(input, 0);

	if(parameters.potential.type==1){
		parameters.currentDomain = 0;
		solveODE(z, parameters, y1);

		parameters.prevDomainY0=y1[0];
		parameters.prevDomainY1=y1[1];
		parameters.prevDomainY2=y1[2];
		parameters.currentDomain = 1;
		solveODE(z, parameters, y2);

		gsl_vector_set(f, 0, y2[0]);
		gsl_vector_set(f, 1, y2[2]-1);
	}
	else if(parameters.potential.type==2){
		parameters.currentDomain = 0;
		solveODE(z, parameters, y1);

		parameters.prevDomainY0=y1[0];
		parameters.prevDomainY1=y1[1];
		parameters.prevDomainY2=y1[2];
		parameters.currentDomain = 1;
		solveODE(z, parameters, y2);

		parameters.prevDomainY0=y2[0];
		parameters.prevDomainY1=y2[1];
		parameters.prevDomainY2=y2[2];
		parameters.currentDomain = 2;
		solveODE(z, parameters, y3);

		parameters.currentDomain = 0;

		gsl_vector_set(f, 0, y3[0]);
		gsl_vector_set(f, 1, y3[2]-1);
	}
	else{
		fprintf(stderr, "ERROR: incorrect parameters in solveODEMultipleDomains()\n");
		exit(EXIT_FAILURE);
	}
	return GSL_SUCCESS;
}

// Use GSL to find (E,z) so that the final conditions are met (phi(L)=0 & N(L)=1)
void findMultipleRoots(schrodingerParameters params, double roots[2]){
	gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (gsl_multiroot_fsolver_hybrids, 2);
	int status;
	size_t i, iter = 0;
	gsl_multiroot_function f = {&solveODEMultipleDomains, 2, &params};
	double x_init[2] = {1.0, 1.0};
	gsl_vector *x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, x_init[0]);
	gsl_vector_set (x, 1, x_init[1]);
	
	gsl_multiroot_fsolver_set (s, &f, x);

	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate(s);

		if (status)
			break;

		status = gsl_multiroot_test_residual (s->f, 1e-8);
	}
	while (status == GSL_CONTINUE && iter < 1000);

	roots[0]=gsl_vector_get (s->x, 0); //E
	roots[1]=gsl_vector_get (s->x, 1); //z

	printf ("status = %s\n", gsl_strerror (status));

	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (x);
}

// Write in a file the points used to draw V(x)
void savePotential(schrodingerParameters params){
	FILE* dataFile = fopen("data/potential.dat", "w");
	fprintf(dataFile, "x V(x)\n");
	double step=0.0;
	if(params.potential.type==0){ // because there is only one value so 2 points are enough to draw the line
		step=params.bound/2;
	}
	else{
		step=params.bound/10;
	}

	double potentialX=0.0;
	double prev_potentialX=0.0;
	for(double x=0.0; x<=params.bound; x+=step){
		potentialX = getPotential(x, params.potential);
		// to get a vertical line and not a leaning wall
		if(prev_potentialX<potentialX){
			fprintf(dataFile, "%lf %lf\n", x, 0.0);
		}
		else if(prev_potentialX>potentialX){
			fprintf(dataFile, "%lf %lf\n", x-step, 0.0);
		}
		fprintf(dataFile, "%lf %lf\n", x, potentialX);
		prev_potentialX=potentialX;
	}
	fclose(dataFile);
}

//Solve Schrodinger's equation with the given parameters
void solveSchrodinger(schrodingerParameters* params){
	if(params->potential.type==0){
		double z;
		z=findRoot(*params);
		params->doDraw=1;
		if(solveODE(z, *params, NULL) == GSL_SUCCESS)
			printf("SUCCESS: the equation was solved\n");
	}
	else if(params->potential.type==1 || params->potential.type==2){
		double z=0.0;
		double y[3]={0.0, 0.0, 0.0};
		double roots[2] = {0.0, 0.0};

		findMultipleRoots(*params, roots);
		params->energy=roots[0];
		params->doDraw=1;
		for(int i_domain=0; i_domain<params->potential.type+1; i_domain++){
			params->currentDomain=i_domain;
			if(solveODE(roots[1], *params, y) != GSL_SUCCESS){
				fprintf(stderr, "ERROR : in solveSchrodinger(), the ODE n°%d could not be solved\n", i_domain);
				exit(EXIT_FAILURE);
			}
			params->prevDomainY0=y[0];
			params->prevDomainY1=y[1];
			params->prevDomainY2=y[2];
		}
		printf("SUCCESS: the equation was solved\n");
	
	}
	else{
		fprintf(stderr, "ERROR : in solveSchrodinger(), parameters are incorrect\n");
		exit(EXIT_FAILURE);
	}
	savePotential(*params);

}