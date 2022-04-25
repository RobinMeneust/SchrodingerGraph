#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

/*
	
alpha = (2m/hbar²)*(E-V(x))
Later on we will have:
params[0] = 2m/hbar²
params[1] = E

y(x) = {phi(x),phi1(x), N(x)}
f(x,y) = {phi1, alpha*phi, phi*phi}			(phi1=(dphi/dx))

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


// used to calculate alpha, which is Schrodinger's aquation : phi'' + alpha*phi = 0
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

int main(){
	double z=1; //temp value
	double constants[2] = {1, 1} ; //(temporary values) params[0] = 2m/hbar²	params[1] = E
	gsl_odeiv2_system sys = {func, NULL, 3, &constants}; // we initialize the ODE system
	gsl_odeiv2_driver * driver = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0.0); 
	double x = 0.0, x1 = 10.0; // we set the bounds
	double y[3] = {0 , z, 0}; //(y0,y1,y2)(t=0) = (0,z,0)
	
	//We get the values of (y0,y1,y2)
	for (int i = 1; i <= 10; i++)
	{
		double xi = i * x1 / 10.0;
		int status = gsl_odeiv2_driver_apply (driver, &x, xi, y);
		if (status != GSL_SUCCESS)
		{
			printf ("error, return value=%d\n", status);
			break;
		}
		printf ("%.5e %.5e %.5e %.5e\n", x, y[0], y[1], y[2]); // We print the values, it will be saved in a file to be used by gnuplot in a future version
	}
	
	gsl_odeiv2_driver_free (driver);
	return 0;
}