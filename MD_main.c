/*
 MD_main.c
 
 Created by AL on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "initfcc.h"
#include "alpotential.h"


/* Main program */
int main()
{

    	/* Declaration of variables */
	double latticeSpacing;
	double timestep;
	double nbr_of_timesteps;
	double nbr_of_atoms;
	double energy;
	double pe;
	double ke;


	/* Initiation of variables */
	latticeSpacing = 1; //Should be in metal units
	timestep = 1;
	nbr_of_timesteps = 1000;
	nbr_of_atoms = 256;


	// Initiate a fcc lattice of Al-particles

	// Introduce displacements, about 5 % of the lattice spacing

	// Set initial velocities to 0 

	// Save the time evolution of the pe, ke and total energy
	
	/* Print energy data to output file */
	FILE *e_file;
	e_file = fopen("energy.data","w");

	fprintf(e_file,"%.5f \t %e \t %e \n", energy, pe, ke);

	fclose(e_file);
    
}
