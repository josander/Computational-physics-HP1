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
	double lattice_param;
	double timestep;
	double nbr_of_timesteps;
	double nbr_of_atoms;
	int Nx, Ny, Nz;
	double energy;
	double pe;
	double ke;
     	double random_value;

	/* Initiation of variables */
	lattice_param = 4.05; // Units: [Å]
	timestep = 0.1;
	nbr_of_timesteps = 1000;
	nbr_of_atoms = 256;
	Nx = 4, Ny = 4, Nz = 4;

	/* Declaration of matrixes */
	double positions[4*Nx*Ny*Nz][3];

	/* Initiation of the fcc lattice of Al-particles */
	init_fcc(positions, nbr_of_atoms, lattice_param);

	// Introduction of displacements, about 5 % of the lattice spacing
     	srand(time(NULL)); // Should only be called once
     	random_value = (double) rand() / (double) RAND_MAX;

	// Set initial velocities to 0 

	// Use the Verlet algorithm

	// Save the time evolution of the pe, ke and total energy
	
	/* Print energy data to output file */
	FILE *e_file;
	e_file = fopen("energy.data","w");

	fprintf(e_file,"%.5f \t %e \t %e \n", energy, pe, ke);

	fclose(e_file);
    
}
