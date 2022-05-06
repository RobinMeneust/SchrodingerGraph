/**
 * \file schrodinger_functions.c
 * \brief Contains functions used to solve the Schrodinger equation with the given parameters
 * \date 2022
 */

#include "../include/types_and_constants.h"
#include "../include/schrodinger_functions.h"

/**
 * \fn double getPotential(double x, potentialParams potential)
 * \brief Get V(x) at a given x and with the given parameters that define V(x)
 * \param x Point where we want to get the potential
 * \param potential Contains all the information about the potential, it gives its value
 * \return Value of the potential at the point x: V(x)
 */

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

/**
 * \fn double getCoefficientFromParams(double x, schrodingerParameters params)
 * \brief Calculate the coefficient of phi which is in the equation of Schrodinger : phi''(x) + C(x)*phi(x) = 0
 * \param x Point where we want to get the coefficient
 * \param params Contains all the parameters sent to the main solver
 * \return Value of the coefficient of phi(x) in the equation
 */

double getCoefficientFromParams(double x, schrodingerParameters params)
{
	return params.alpha*(params.energy-getPotential(x, params.potential));
}

/**
 * \fn int y_derivative (double x, const double y[], double f[], void *params)
 * \brief Calculate f(y0,y1,y2) which is the derivative of the vector (y0(x),y1(x),y2(x))
 * \param x Point where we want to get the derivative
 * \param y Previous value of the vector (y0(x),y1(x),y2(x))
 * \param f Derivative f(y0,y1,y2)
 * \param params Contains all the parameters sent to the main solver
 * \return Status: it's GSL_SUCCESS if there was no error
 */

int y_derivative (double x, const double y[], double f[], void *params){
	double coeff = getCoefficientFromParams(x, *(schrodingerParameters*) params);
	f[0]=y[1];
	f[1]=(-1.0)*coeff*y[0];
	f[2]=y[0]*y[0];

	return GSL_SUCCESS;
}

/**
 * \fn int solveODE(double z, schrodingerParameters params, double f[3])
 * \brief Solve the ODE with the given parameters to get the function phi(x), and draw it if params.doDraw!=0
 * \param z Derivative of phi(x) at x=0 (i.e: z = y1(0))
 * \param params Contains all the parameters sent to the main solver
 * \param f Value of the vector (y0(L),y1(L),y2(L)) (L is the last value taken by x)
 * \return Status: it's GSL_SUCCESS if there was no error
 */

int solveODE(double z, schrodingerParameters params, double f[3]){
	FILE* dataFile=NULL;

	double x=0.0;
	double length=params.bound;
	double y1_ini=0.0, y2_ini=z, y3_ini=0.0; // phi(0)=0 and N(0)=0

	if(params.doDraw){
		dataFile = fopen("data/phi.dat", "a");
		if(dataFile==NULL){
			fprintf(stderr, "ERROR: the file data/phi.dat can't be opened or created. Please check if the folder data exists\n");
			exit(EXIT_FAILURE);
    		}
	}
	
	gsl_odeiv2_system sys = {y_derivative, NULL, 3, &params}; // we initialize the ODE system
	// we set the bounds

	if(params.currentDomain==0)
	{
		if(params.potential.type!=0)
			length=params.potential.a-x;
	}
	else{
		y1_ini=params.prevDomainY0;
		y2_ini=params.prevDomainY1;
		y3_ini=params.prevDomainY2;
		if(params.currentDomain==1){
			x = params.potential.a;
			length=params.potential.b-x;
		}
		else if(params.currentDomain==2){
			x = params.potential.b;
			length=params.bound-x;
		}
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
		int status = gsl_odeiv2_driver_apply(driver, &x, xi, y);
		if (status != GSL_SUCCESS)
		{
			printf ("ERROR: in solveODE() %d\n", status);
			break;
		}
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

/**
 * \fn double findRoot(schrodingerParameters params)
 * \brief Find the value of z=y1(0) so that the edge conditions are met
 * \param params Contains all the parameters sent to the main solver
 * \return y1(0)
 */

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

/**
 * \fn int solveODEMultipleDomains(const gsl_vector* input, void* params, gsl_vector* f)
 * \brief Check if (E,z) allow us to meet the edge conditions. It Calls solveODE() to solve the equation on different domains with different potentials
 * \param input Vector (E,z) that we are cheking
 * \param params Contains all the parameters sent to the main solver
 * \param f Value of the vector (y0(L),y1(L),y2(L)) (L is the last value taken by x)
 * \return Status: it's GSL_SUCCESS if there was no error
 */

int solveODEMultipleDomains(const gsl_vector* input, void* params, gsl_vector* f){
	schrodingerParameters parameters = *(schrodingerParameters*) params;
	double y[3]={0.0};
	double z=gsl_vector_get(input, 1);
	parameters.energy=gsl_vector_get(input, 0);

	if(parameters.energy<=0){ // a negative energy caused the program to freeze
		gsl_vector_set(f, 0, 1.0);
		gsl_vector_set(f, 1, 1.0);
		return GSL_SUCCESS;
	}

	if(parameters.potential.type==1 || parameters.potential.type==2){
		for(int i_domain=0; i_domain<parameters.potential.type+1; i_domain++){
			parameters.currentDomain = i_domain;
			if(solveODE(z, parameters, y) != GSL_SUCCESS){
				fprintf(stderr, "ERROR : in solveODEMultipleDomains(), the ODE n°%d could not be solved\n", i_domain);
				exit(EXIT_FAILURE);
			}
			parameters.prevDomainY0=y[0];
			parameters.prevDomainY1=y[1];
			parameters.prevDomainY2=y[2];
		}
		gsl_vector_set(f, 0, y[0]);
		gsl_vector_set(f, 1, y[2]-1);
	}
	else{
		fprintf(stderr, "ERROR: incorrect parameters in solveODEMultipleDomains()\n");
		exit(EXIT_FAILURE);
	}
	return GSL_SUCCESS;
}

/**
 * \fn int findMultipleRoots(double x_init, schrodingerParameters params, double roots[2])
 * \brief Find the roots (E,z) that meet the edge conditions (phi(L)=0 & N(L)=1).
 * \param x_init Initial root value
 * \param params Contains all the parameters sent to the main solver
 * \param roots Vector (E,z) that meet the edge conditions, if it exists
 * \return Status: it's GSL_SUCCESS if there was no error
 */

int findMultipleRoots(double x_init, schrodingerParameters params, double roots[2]){
	gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (gsl_multiroot_fsolver_hybrid, 2);
	int status=GSL_CONTINUE;
	size_t iter = 0;
	gsl_multiroot_function f = {&solveODEMultipleDomains, 2, &params};
	gsl_vector *x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, x_init);
	gsl_vector_set (x, 1, x_init);
	
	gsl_multiroot_fsolver_set (s, &f, x);
	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate(s);
		if (status)
			break;
		status = gsl_multiroot_test_residual (s->f, 1e-8);

	}while (status == GSL_CONTINUE && iter < 1000);

	roots[0]=gsl_vector_get (s->x, 0); //E
	roots[1]=gsl_vector_get (s->x, 1); //z

	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (x);

	return status;
}

/**
 * \fn void savePotential(schrodingerParameters params)
 * \brief Write in a file the points used to draw V(x)
 * \param params Contains all the parameters sent to the main solver
 */

void savePotential(schrodingerParameters params){
	FILE* dataFile = fopen("data/potential.dat", "w");
	if(dataFile==NULL){
		fprintf(stderr, "ERROR: the file data/potential.dat can't be opened or created. Please check if the folder data exists\n");
		exit(EXIT_FAILURE);
    	}
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

/**
 * \fn void solveSchrodinger(schrodingerParameters* params)
 * \brief Solve Schrodinger's equation with the given parameters
 * \param params Contains all the parameters used by this function to solve the equation
 */

void solveSchrodinger(schrodingerParameters* params){
	double z=0.0;
	if(params->potential.type==0){
		z=findRoot(*params);
		params->doDraw=1;
		if(solveODE(z, *params, NULL) == GSL_SUCCESS)
			printf("SUCCESS: the equation was solved\n");
	}
	else if(params->potential.type==1 || params->potential.type==2){
		double y[3]={0.0, 0.0, 0.0};
		double roots[10][2]={{0.0}};
		int isAlreadyInArray;
		int i_min_energy=0.0;
		int i=0;
		// we get the possible values for (E,z)
		for(double x=0; x<5; x+=0.1){
			if(findMultipleRoots(x, *params, roots[i])==GSL_SUCCESS){
				isAlreadyInArray=0;
				// We don't want to have multiple occurrences of one value so we don't save it in that case
				for(int j=0; j<i; j++){
					if(fabs(roots[i][0]-roots[j][0])<1e-5){ // if they are too close we ignore it
						isAlreadyInArray=1;
						break;
					}
				}
				if(!isAlreadyInArray){
					i++;
					if(i>=10)
						break;
				}
			}
		}
		if(i==0){ // if we haven't found any solution
			fprintf(stderr, "ERROR : in solveSchrodinger(), no solution was found\n");
			exit(EXIT_FAILURE);
		}
		// we search the minimum energy level Ei (here we can choose any energy level, it's just an arbitrary choice)
		i_min_energy = 0;
		for(i=1; i<10; i++){
			if(roots[i][0]!=0 && roots[i][0] < roots[i_min_energy][0])
				i_min_energy=i;
		}

		params->energy=roots[i_min_energy][0];
		params->doDraw=1;

		// We solve the equation in the 2 or 3 domains
		for(int i_domain=0; i_domain<params->potential.type+1; i_domain++){
			params->currentDomain=i_domain;
			if(solveODE(roots[i_min_energy][1], *params, y) != GSL_SUCCESS){
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
