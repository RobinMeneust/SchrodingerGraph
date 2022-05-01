#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>

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
	double prevDomainY2; // N(x) at the right bound of the previous domain
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
int func (double x, const double y[], double f[] /* = dydt*/, void *params){
	double alpha = getCoefficientFromParams(x, *(schrodingerParameters*) params);
	f[0]=y[1];
	f[1]=(-1.0)*alpha*y[0];
	f[2]=y[0]*y[0];

	return GSL_SUCCESS;
}

// solve the ODE with the giveen parameters to get the function phi(x), and draw it if params.doDraw!=0

int solveODE(double z, schrodingerParameters params, double f[2]){
	FILE* dataFile=NULL;
	double x=0.0;
	double length=params.bound;
	double y1_ini=0.0, y3_ini=0.0; // phi(0)=0 and N(0)=0

	if(params.doDraw)
		dataFile = fopen("data/phi.dat", "a");
	
	gsl_odeiv2_system sys = {func, NULL, 3, &params}; // we initialize the ODE system
	// we set the bounds
	if(params.potential.type!=0 && params.currentDomain==0){
		x = 0.0;
		length=params.potential.a-x;
	}
	else if(params.currentDomain==1){
		x = params.potential.a;
		length=params.potential.b-x;
		y1_ini=params.prevDomainY0;
		y3_ini=params.prevDomainY2;
	}
	else if(params.currentDomain==2){
		x = params.potential.b;
		length=params.bound-x;
		y1_ini=params.prevDomainY0;
		y3_ini=params.prevDomainY2;
	}

	gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0);

	double y[3] = {y1_ini, z, y3_ini}; //(y0,y1,y2)(x_ini) = (0,z,0)
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
		f[0]=y[0]; //y0(L)
		f[1]=y[2]-1; //y3(L)-1
	}

	gsl_odeiv2_driver_free (driver);
	if(params.doDraw)
		fclose(dataFile);
	
	return GSL_SUCCESS;
}

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
		solveODE(z, params, y);
		//if(fabs(y[0])<0.1 && fabs(y[1])<0.1)
		//	printf("y0(L)= %lf | y3(L)-1 = %lf\n", y[0], y[1]);
		eps=fabs(y[1]);
		if(eps<eps_min){
			eps_min=eps;
			best_z_approx=z;
		}
		iter++;
	}while(eps>1e-8 && iter<5000);
	
	return best_z_approx;
}

void findMultipleRoots(schrodingerParameters params, double roots[3]){
	schrodingerParameters paramsSubdivided=params; // params for one of the 3 domains ([0,a] [a,b] [b,L])
	double y1[2]={0.0, 0.0}; // Y(a) in 1st domain
	double y2[2]={0.0, 0.0}; // Y(b) in 2nd domain
	double y3[2]={0.0, 0.0}; // Y(L) in 3rd domain

	double step=0.01;
	double z[3]={0.0, 0.0, 0.0};
	double eps=3.0;
	int iter=0;
	int max_iter=10000;
	double eps_min=eps;
	if(params.potential.type==2){
		do{
			z[0]+=step;
			do{
				z[1]+=step;
				do{
					z[2]+=step;
					//printf("Z %lf %lf %lf\n", z[0], z[1], z[2]);
					params.currentDomain = 0;
					solveODE(z[0], params, y1);

					params.prevDomainY0=y1[0];
					params.prevDomainY2=y1[2];
					params.currentDomain = 1;
					solveODE(z[1], params, y2);

					params.prevDomainY0=y2[0];
					params.prevDomainY2=y2[2];
					params.currentDomain = 2;
					solveODE(z[2], params, y3);

					params.currentDomain = 0;
					//if(fabs(y[0])<0.1 && fabs(y[1])<0.1)
					//	printf("y0(L)= %lf | y3(L)-1 = %lf\n", y[0], y[1]);
					eps=fabs(y1[1])+fabs(y2[1])+fabs(y3[1]);
					if(eps<eps_min){
						//printf("iter %d | Z=%lf %lf %lf | EPS %lf\n", iter, z[0], z[1], z[2], eps);
						eps_min=eps;
						roots[0]=z[0];
						roots[1]=z[1];
						roots[2]=z[2];
					}
					iter++;
				}while(z[2]<5.0 && iter<max_iter);
				z[2]=0;
			}while(z[1]<5.0  && iter<max_iter);
			z[1]=0;
		}while(eps>1e-8 && iter<max_iter);
		printf("iter %d | Z=%lf %lf %lf | EPS %lf\n", iter, z[0], z[1], z[2], eps);
	}
	else{
		fprintf(stderr, "ERROR: not yet implemented\n");
		exit(EXIT_FAILURE);
	}
	printf("ROOTS %lf %lf %lf | ACCURACY %lf\n", roots[0], roots[1], roots[2], eps_min);
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
	if(params.potential.type==0){
		double z;
		z=findRoot(params);
		params.doDraw=1;
		if(solveODE(z, params, NULL) == GSL_SUCCESS)
			printf("SUCCESS: the equation was solved\n");
	}
	else if(params.potential.type==1 || params.potential.type==2){
		double z[3]={0.0, 0.0, 0.0};
		double y[3]={0.0, 0.0, 0.0};
		findMultipleRoots(params, z);
		params.doDraw=1;
		for(int i_domain=0; i_domain<params.potential.type+1; i_domain++){
			params.currentDomain=i_domain;
			if(solveODE(z[i_domain], params, y) != GSL_SUCCESS){
				fprintf(stderr, "ERROR : in solveSchrodinger(), the ODE n°%d could not be solved\n", i_domain);
				exit(EXIT_FAILURE);
			}
			params.prevDomainY0=y[0];
			params.prevDomainY2=y[2];
		}
		printf("SUCCESS: the equation was solved\n");
			
	}
	else{
		fprintf(stderr, "ERROR : in solveSchrodinger(), parameters are incorrect\n");
		exit(EXIT_FAILURE);
	}

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

	FILE* dataFile = fopen("data/phi.dat", "w"); fprintf(dataFile, "x phi(x)\n"); fclose(dataFile); // We initialize the file phi.dat

	double m = 6.25;
	double l = 1.0;
	double hbar = 0.65625;
	double energy = hbar*hbar*4*M_PI*M_PI / (8*m*l*l); // ~=0.34 according to the formula : En = h^2 * n^2 / (8*m*l^2) 		So to get the n-th energy level we can just multiply by n^2 this value
	double alpha = 2*m / (hbar*hbar); // we have (d^2)(phi(x))/d(x^2) + alpha phi(x) = 0
	potentialParams potential;
	schrodingerParameters params;

	params.energy=energy;
	params.bound=l;
	params.potential.a=0.0;
	params.potential.b=l;
	params.alpha=alpha;
	params.doDraw=0;
	params.currentDomain=0;
	params.prevDomainY0=0.0;
	params.prevDomainY2=0.0;

	potential.type=0;
	potential.v0=1.0;
	
	printf("What case do you want to view\n0: V(x)=0 everywhere\n1: V(x) is a step\n2: V(x) is rectangular\nANSWER: "); scanf("%d", &(potential.type));
	if(potential.type!=0){ // if V(x) != 0 for all x
		printf("Give the value of the potential in the well when it's not null: V0 = "); scanf("%lf", &(potential.v0));
		potential.a = 0.4;
		potential.b = 0.6;
	}

	params.potential=potential;
	
	
	//FILE* dataFile = fopen("data/phiALL.dat", "w");fclose(dataFile);
	solveSchrodinger(params);


	return 0;
}