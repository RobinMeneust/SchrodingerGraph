#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multiroots.h>

#define N_POINTS 100 //number of points
#define L 10.0

typedef struct potentialParams{
	int type;
	double value;
}potentialParams;

/*
	
alpha = (2m/hbar²)*(E-V(x))
Later on we will have:
params[0] = 2m/hbar²
params[1] = E

y(x) = {phi(x),phi1(x), N(x)}
f(x,y) = {phi1, alpha*phi, phi*phi}			(phi1=(dphi/dx))

y(x)=0 <=> x=?

d(phi)/dx - alpha * phi = 0

*/


// It will change between the 3 cases that we will see here (rectangular, null...)
double getPotential(double x){
	return 0;
	/*
	if(x>5){
		...
	}
	else{
		...
	}
	*/
}


// Calculate alpha which is in the equation of Schrodinger : phi'' + alpha*phi = 0
double getConstantFromParams(double x, double* params)
{
	return params[0]*(params[1]-getPotential(x));
}

// calculate f(y0,y1,y2)
int func (double x, const double y[], double f[] /* = dydt*/, void *params){
	double alpha = getConstantFromParams(x, (double*) params);
	f[0]=y[1];
	f[1]=(-1.0)*alpha*y[0];
	f[2]=y[0]*y[0];

	return GSL_SUCCESS;
}

// get an approximation of the values of the function phi(x) in the equation, to find the root
int functionForRoot(const gsl_vector * input, void *params, gsl_vector * f)
{
	const double energy = gsl_vector_get (input, 0);
	const double z = gsl_vector_get (input, 1);

	double constants[2] = {1, energy} ; //(temporary values) params[0] = 2m/hbar²	params[1] = E
	gsl_odeiv2_system sys = {func, NULL, 3, &constants}; // we initialize the ODE system
	double x = 0.0, l = 10.0; // we set the bounds

	gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0); 

	double y[3] = {0, z, 0}; //(y0,y1,y2)(t=0) = (0,z,0) if V(x)=0 everywhere for example
	
	//We get the values of (y0,y1,y2)
	FILE* dataFile = fopen("data/phiALL.dat", "a"); // USED FOR TESTING
	for (int i = 1; i <= N_POINTS; i++)
	{
		double xi = i * l / (double)N_POINTS;
		int status = gsl_odeiv2_driver_apply (driver, &x, xi, y);
		if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}
		if(y[0]<100 && y[0]>-100 && x<100 && x>-100){
			fprintf(dataFile, "%lf %lf\n", x, y[0]);
		}
	}
	fclose(dataFile);
	gsl_vector_set (f, 0, y[0]); //y0(L)
	gsl_vector_set (f, 1, y[2]-1); //y2(L)-1

	gsl_odeiv2_driver_free (driver);

	return GSL_SUCCESS;
}

//FOR DEBUG
/*
int print_state (size_t iter, gsl_multiroot_fsolver * s)
{
	printf ("iter = %3lu x = % .3f % .3f | f(x) = % .3e % .3e\n", 
		iter,
		gsl_vector_get (s->x, 0),
		gsl_vector_get (s->x, 1),
		gsl_vector_get (s->f, 0),
		gsl_vector_get (s->f, 1)
	);
}
*/
// Find the roots for the shooting method
void findRoots(double* roots){
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;

	int status;
	size_t i, iter = 0;

	const size_t n = 2;
	gsl_multiroot_function f = {&functionForRoot, n, NULL};

	double x_init[2] = {-5.0, -5.0};
	gsl_vector *x = gsl_vector_alloc (n);

	gsl_vector_set (x, 0, x_init[0]);
	gsl_vector_set (x, 1, x_init[1]);

	T = gsl_multiroot_fsolver_hybrids;
	s = gsl_multiroot_fsolver_alloc (T, 2);
	gsl_multiroot_fsolver_set (s, &f, x);

	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate(s);

		if (status)
			break;

		status = gsl_multiroot_test_residual (s->f, 1e-7);
	}
	while (status == GSL_CONTINUE && iter < 1000);

	roots[0]=gsl_vector_get (s->x, 0); //E
	roots[1]=gsl_vector_get (s->x, 1); //z

	printf ("status = %s\n", gsl_strerror (status));

	gsl_multiroot_fsolver_free (s);
	gsl_vector_free (x);
}


// solve the ODE with the correct params (that we got from findRoots())
void solveAndSaveData(double* params){
	FILE* dataFile = fopen("data/phi.dat", "w");
	
	const double energy = params[0];
	const double z = params[1];

	double constants[2] = {1, energy} ; //(temporary values) params[0] = 2m/hbar²	params[1] = E
	gsl_odeiv2_system sys = {func, NULL, 3, &constants}; // we initialize the ODE system
	double x = 0.0, l = 10.0; // we set the bounds

	gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0); 

	double y[3] = {0 , z, 0}; //(y0,y1,y2)(t=0) = (0,z,0)
	fprintf(dataFile, "x phi(x)\n");
	fprintf(dataFile, "%lf %lf\n", x, y[0]);
	//We get the values of (y0,y1,y2)
	for (int i = 1; i <= N_POINTS; i++)
	{
		double xi = i * l / (double)N_POINTS;
		int status = gsl_odeiv2_driver_apply (driver, &x, xi, y);
		if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}
		//printf("f(%.2lf)=%lf\n", x, y[0]);
		fprintf(dataFile, "%lf %lf\n", x, y[0]);
	}

	gsl_odeiv2_driver_free (driver);

	fclose(dataFile);
}

void drawPotential(){
	FILE* dataFile = fopen("data/potential.dat", "w");
	fprintf(dataFile, "x V(x)\n");
	for(double x=0; x<=L; x++){
		fprintf(dataFile, "%lf %lf\n", x, getPotential(x));
	}
	fclose(dataFile);
}

int main(){
	double roots[2]={0,0};
	FILE* dataFile = fopen("data/phiALL.dat", "w");fclose(dataFile);
	/*
	potentialParams potential;
	potential.type=0;
	potential.value=1; //temporary value
	
	printf("What case do you want to view\n0: V(x)=0 everywhere\n1: V(x) is a step\n2: V(x) is rectangular\n: "); scanf("%d", &(potential.type));
	*/
	findRoots(roots);
	solveAndSaveData(roots);
	drawPotential();

	return 0;
}