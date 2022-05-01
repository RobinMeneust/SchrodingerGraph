#include <stdio.h>
#include <stdlib.h>
#include "../include/schrodinger_functions.h"
#include "../include/types_and_constants.h"

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
	params.prevDomainY1=0.0;
	params.prevDomainY2=0.0;

	potential.type=0;
	potential.v0=1.0;
	
	printf("What case do you want to view\n0: V(x)=0 everywhere\n1: V(x) is a step\n2: V(x) is rectangular\nANSWER: "); scanf("%d", &(potential.type));

	if(potential.type<0 || potential.type>2){
		fprintf(stderr, "ERROR: incorrect value. You must provide a number between 0 and 2\n");
		exit(EXIT_FAILURE);
	}

	if(potential.type!=0){ // if V(x) != 0 for all x
		printf("Give the value of the potential in the well when it's not null: V0 = "); scanf("%lf", &(potential.v0));
		potential.a = 0.4;
		if(potential.type==1)
			potential.b = l; // because we only have 2 domains: [0,a] and [a,l]
		else
			potential.b=0.6;
	}

	params.potential=potential;
	
	
	solveSchrodinger(params);


	return 0;
}