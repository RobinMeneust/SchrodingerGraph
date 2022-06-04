/**
 * \file main.c
 * \brief It gets the user's input, checks it, initializes the variables and sends it to the solver to get the required data to draw the graph of a solution to the Schrodinger equation
 * \date 2022
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/schrodinger_functions.h"
#include "../include/types_and_constants.h"

#define MASS_PROTON 10456.25;
#define MASS_ELECTRON 5.693;

/*
We choose those units to get values near to 1 :
length = 1e-9 m = 1 nm
energy :  eV
mass:
	energy: 1 J = (1/1.6) * 1e19 eV
	length: 1 m = 1e9 nm => m^-2 = 1e-18 nm^-2
	time: 1 s = 1e15 fs => (1s)^2 = 1e30 fs^2		(fs = "femtosecond")

	and the mass (kg) is homogeneous to E/c^2 (J.s^2/m^2) so we have:

	1 kg 	= 1 J.s².m^-2
			= (1/1.6) * 1e(19-18) eV.s^2.nm^-2
			= (1/1.6) * 10 eV.s^2.nm^-2
			= (1/1.6) * 10 * 1e30 eV.fs^2.nm^-2
			= 6.25 * 1e30 eV.nm^-2.fs^2

	and m  	= 9.109 * 1e-31 kg
			= 9.109 * 1e-31 * 6.25 * 10^30 eV.nm^-2.fs^2
			= 9.109 * 0.1 * 6.25 eV.nm^-2.fs^2
			= 5.693 eV.nm^-2.fs^2


hbar 	= 1.05e-34 J.s 
		= 1.05e-34 * (1/1.6) * 1e19 * 1e15 eV.fs
		= 0.65625 * 1e(-34+19+15) eV.fs
		= 0.65625 eV.fs
*/

/**
 * \fn int main()
 * \brief Main function of this project
 * \return 0 if the function runs and exits correctly
 */

int main(){
	// We initialize the file phi.dat
	FILE* dataFile = fopen("data/phi.dat", "w"); 
	FILE* dataFileSquare = fopen("data/phi_square.dat", "w"); 
	if(dataFile==NULL || dataFileSquare==NULL){
        fprintf(stderr, "ERROR: the file data/phi.dat can't be opened or created. Please check if the folder data exists\n");
        exit(EXIT_FAILURE);
    }
	fprintf(dataFile, "x phi(x)\n");fclose(dataFile);
	fprintf(dataFileSquare, "x |phi(x)|²\n");fclose(dataFileSquare);

	double m = MASS_ELECTRON;
	double l = 1.0;
	double hbar = 0.65625;
	double energy = hbar*hbar*4*M_PI*M_PI / (8*m*l*l); // ~=0.34 according to the formula : En = h^2 * n^2 / (8*m*l^2) 		So to get the n-th energy level we can just multiply by n^2 this value
	double alpha = 2*m / (hbar*hbar); // we have (d^2)(phi(x))/d(x^2) + alpha (E-V) phi(x) = 0
	int energyLevel=1;
	char potentialTypeInText[50];
	potentialParams potential;
	schrodingerParameters params;

	params.energy=energy;
	params.bound=l;
	params.alpha=alpha;
	params.doDraw=0;
	params.currentDomain=0;
	params.prevDomainY0=0.0;
	params.prevDomainY1=0.0;
	params.prevDomainY2=0.0;

	potential.type=0;
	potential.v0=1.0;
	potential.a=0.0;
	potential.b=l;
	
	printf("What case do you want to view\n0: V(x)=0 everywhere\n1: V(x) is a step\n2: V(x) is rectangular\n3: V(x) is parabolic\nANSWER: "); scanf("%d", &(potential.type));

	if(potential.type<0 || potential.type>3){
		fprintf(stderr, "ERROR: incorrect value. You must provide a number between 0 and 2\n");
		exit(EXIT_FAILURE);
	}

	if(potential.type!=0){ // if V(x) != 0 for all x
		printf("Give the value of the potential in the well when it's not null, between 0 and 60: V0 = "); scanf("%lf", &(potential.v0));
		if(potential.v0<0 || potential.v0>60){
			fprintf(stderr, "ERROR: You need to provide a value in [0,60] for V0\n");
			exit(EXIT_FAILURE);
		}
		potential.a = 0.45;
		if(potential.type==1)
			potential.b = l; // because we only have 2 domains: [0,a] and [a,l]
		else
			potential.b = 0.55;
	}
	else{
		printf("Give the energy level that you want (n>0 and n is an integer): n = "); scanf("%d", &energyLevel);
		if(energyLevel>0)
			params.energy*= (energyLevel*energyLevel);
		else{
			fprintf(stderr, "ERROR: You need to provide a strictly positive integer\n");
			exit(EXIT_FAILURE);
		}
	}
	//Used by gnuplot
	switch (potential.type)
	{
		case 0:
			strcpy(potentialTypeInText, "null potential in an infinite well");
			break;
		case 1:
			strcpy(potentialTypeInText, "step potential in an infinite well");
			break;
		case 2:
			strcpy(potentialTypeInText, "rectangular potential in an infinite well");
			break;
		case 3:
			strcpy(potentialTypeInText, "parabolic potential in an infinite well");
			break;
		default:
			fprintf(stderr, "ERROR: incorrect potential type in main()\n");
			exit(EXIT_FAILURE);
	}

	params.potential=potential;
	
	solveSchrodinger(&params);

	dataFile = fopen("data/phi.dat", "a");
	if(dataFile==NULL){
        fprintf(stderr, "ERROR: the file data/phi.dat can't be opened or created. Please check if the folder data exists\n");
        exit(EXIT_FAILURE);
    }
	if(params.energy<0.01)
		fprintf(dataFile, "\npotential_type & energy:\n%s\n%.3e", potentialTypeInText, params.energy);
	else
		fprintf(dataFile, "\npotential_type & energy:\n%s\n%.3lf", potentialTypeInText, params.energy);
	fclose(dataFile);

	return 0;
}
