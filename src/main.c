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
	double a;
	double b;
}potentialParams;

typedef struct schrodingerParameters{
	potentialParams potential;
	double energy;
	double bound;
	double alpha;
}schrodingerParameters;

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
			fprintf(stderr, "ERROR : in getPotential(), parameters are incorrect");
			exit(EXIT_FAILURE);
	}

	/*
	
	*/
}


// Calculate the coefficient of phi (C) which is in the equation of Schrodinger : phi'' + C*phi = 0
double getCoefficientFromParams(double x, schrodingerParameters params)
{
	return params.alpha*(params.energy-getPotential(x, params.potential));
}

// calculate f(y0,y1,y2)
int func (double x, const double y[], double f[] /* = dydt*/, void *params){
	double alpha = getCoefficientFromParams(x, *(schrodingerParameters*) params);
	f[0]=y[1];
	f[1]=(-1.0)*alpha*y[0];
	f[2]=y[0]*y[0];

	return GSL_SUCCESS;
}

// get an approximation of the values of the function phi(x) in the equation, to find the root
int functionForRoot(double z, schrodingerParameters params, double f[2])
{
	gsl_odeiv2_system sys = {func, NULL, 3, &params}; // we initialize the ODE system
	double x = 0.0, l = params.bound; // we set the bounds

	gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0); 

	double y[3] = {0, z, 0}; //(y0,y1,y2)(t=0) = (0,z,0) if V(x)=0 everywhere for example
	
	//We get the values of (y0,y1,y2)
	//FILE* dataFile = fopen("data/phiALL.dat", "a"); // USED FOR TESTING
	//fprintf(dataFile, "%lf %lf\n", x, y[0]);
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
			//fprintf(dataFile, "%lf %lf\n", x, y[0]);
		}
	}
	//fclose(dataFile);
	f[0]=y[0]; //y0(L)
	f[1]=y[2]-1; //y3(L)-1

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
// Find the root with the shooting method
double findRoot(schrodingerParameters params){
	double y[2]={0.0, 0.0};
	double step=0.01;
	double z=0.0-step;
	double eps=1.0;
	int iter=0;
	double best_z_approx=z;
	double eps_min=eps;

	do{
		z+=step;
		functionForRoot(z, params, y);
		//if(fabs(y[0])<0.1 && fabs(y[1])<0.1)
		//	printf("y0(L)= %lf | y3(L)-1 = %lf\n", y[0], y[1]);
		eps=fabs(y[1]);
		if(eps<eps_min){
			eps_min=eps;
			best_z_approx=z;
		}
		iter++;
	}while(eps>1e-8 && iter<1000);
	
	return best_z_approx;
}


// solve the ODE with the correct params (that we got from findRoots())
void solveAndSaveData(double z, schrodingerParameters params){
	FILE* dataFile = fopen("data/phi.dat", "w");
	
	gsl_odeiv2_system sys = {func, NULL, 3, &params}; // we initialize the ODE system
	double x = 0.0, l = params.bound; // we set the bounds

	gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

	double y[3] = {0, z, 0}; //(y0,y1,y2)(t=0) = (0,z,0)
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

void drawPotential(schrodingerParameters params){
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
		printf("V(%lf) = %lf | prev = %lf \n", x, potentialX, prev_potentialX);
		// to get a vertical line and not a leaning wall
		if(prev_potentialX<potentialX){
			fprintf(dataFile, "%lf %lf\n", x, 0.0);
		}
		else if(prev_potentialX>potentialX){
			fprintf(dataFile, "%lf %lf\n", x, prev_potentialX);
		}
		fprintf(dataFile, "%lf %lf\n", x, potentialX);
		prev_potentialX=potentialX;
	}
	fclose(dataFile);
}

void solveSchrodinger(schrodingerParameters params){
	double z=0.0;
	
	z=findRoot(params);
	solveAndSaveData(z, params);
	drawPotential(params);
}

int main(){
	/*
	We choose those units to get values near to 1 :
	length = 1e-9 m = 1 nm
	energy = 1.6e-19 J = 1 eV
	mass:
		energy: 1 J = (1/1.6) * 1e19 eV
		length: 1 m = 1e9 nm => m^-2 = 1e-18 nm^-2
		time: 1 s = 1e15 fs => (1s)^2 = 1e30 fs^2		(fs = "femtosecond")

		=> J.m^-2.s^2 = (1/1.6) * 10 eV.nm^-2.s^2

		and a mass (kg) is homogeneous to E/c² (J.s²/m²) so we have:
		9,109 × 10−31
		=> m	= 1e-30 kg	= 1e-30 * (1/1.6) * 10 eV.nm^-2.s^2
				= (1/1.6) * 1e-29 eV.nm^-2.s^2
				= (1/1.6) * 10 eV.nm^-2.fs^2
				= 6.25 eV.nm^-2.fs^2

	hbar 	= 1.05e-34 J.s 
			= 1.05e-34 * (1/1.6) * 1e19 * 1e15 eV.fs
			= 0,65625 * 1e(-34+19+15)
			= 0,65625
	*/
	double m = 6.25;
	double l = 1.0;
	double hbar = 0.65625;
	double energy = hbar*hbar*4*M_PI*M_PI / (8*m*l*l); // ~=0.34 according to the formula : En = h^2 * n^2 / (8*m*l^2)
	double alpha = 2*m / (hbar*hbar); // we have (d^2)(phi(x))/d(x^2) + alpha phi(x) = 0

	potentialParams potential;
	potential.type=0;
	potential.v0=1.0;
	
	printf("What case do you want to view\n0: V(x)=0 everywhere\n1: V(x) is a step\n2: V(x) is rectangular\nANSWER: "); scanf("%d", &(potential.type));
	if(potential.type!=0){ // if V(x) != 0 for all x
		printf("Give the value of the potential in the well when it's not null: V0 = "); scanf("%lf", &(potential.v0));
		potential.a = 0.4;
		potential.b = 0.6;
	}

	schrodingerParameters params;
	params.potential=potential;
	params.energy=energy;
	params.bound=l;
	params.alpha=alpha;
	
	//FILE* dataFile = fopen("data/phiALL.dat", "w");fclose(dataFile);
	solveSchrodinger(params);


	return 0;
}